module CtsemCompute
using Distributed, LinearAlgebra, ForwardDiff, Optim, Random, SHA, Serialization, Sockets
using Main.CtsemJulia
const SOURCE_HASH=bytes2hex(sha256(read(@__FILE__)))
const STATE=Ref{Any}(nothing)
getv(x,k)=CtsemJulia.field(x,Symbol(k))
function setup(sd,ids)
    STATE[]=(sd=sd,d=CtsemJulia.prepare_subject_data(sd),ids=collect(ids))
    nothing
end
function blockloss(x,d,ids)
    raw=x[1:d.nparams];tip=x[d.nparams+1:end]
    value=zero(eltype(x))
    for si in ids
        z=copy(raw)
        for j in 1:d.ntipred,p in 1:d.nparams
            ix=d.tie_setup[p,j]
            ix>0 && (z[p]+=tip[ix]*d.tipreds[si,j])
        end
        value-=CtsemJulia.subject_loglik(z,d,d.firstrow[si],d.lastrow[si])
    end
    value
end
function gradientblock(x,d,ids,chunksize)
    fun=z->blockloss(z,d,ids)
    result=ForwardDiff.DiffResults.GradientResult(x)
    cfg=ForwardDiff.GradientConfig(fun,x,ForwardDiff.Chunk{min(chunksize,length(x))}())
    ForwardDiff.gradient!(result,fun,x,cfg)
    (value=ForwardDiff.DiffResults.value(result),gradient=copy(ForwardDiff.DiffResults.gradient(result)))
end
function evaluate(x,chunksize)
    s=STATE[];ids=s.ids;k=min(Threads.nthreads(),length(ids))
    k<=1 && return(gradientblock(x,s.d,ids,chunksize))
    parts=Vector{Any}(undef,k)
    Threads.@threads :static for j in 1:k
        lo=fld((j-1)*length(ids),k)+1;hi=fld(j*length(ids),k)
        parts[j]=gradientblock(x,s.d,ids[lo:hi],chunksize)
    end
    (value=sum(p.value for p in parts),gradient=reduce(+,p.gradient for p in parts))
end
function connect(config)
    ids=Int[]
    try
        for c in config
            n=Int(getv(c,"workers"));n==0 && continue
            flags=["--startup-file=no","--threads="*string(Int(getv(c,"threads"))),"--project="*String(getv(c,"project"))]
            host=String(getv(c,"host"))
            new=if host=="local"
                addprocs(n;exeflags=flags,topology=:master_worker)
            else
                sshflags=Cmd(["-S",String(getv(c,"socket")),"-p",string(Int(getv(c,"port"))),
                    "-o","BatchMode=yes","-o","ConnectTimeout=10"])
                env=["OPENBLAS_NUM_THREADS"=>"1","JULIA_NUM_PRECOMPILE_TASKS"=>"1"]
                depot=String(getv(c,"depot"));isempty(depot)||push!(env,"JULIA_DEPOT_PATH"=>depot)
                addprocs([(host,n)];tunnel=true,topology=:master_worker,sshflags=sshflags,
                    exename=String(getv(c,"julia")),dir=String(getv(c,"source")),exeflags=flags,env=env)
            end
            append!(ids,new)
            path=String(getv(c,"source"));hashes=String.(collect(getv(c,"hashes")))
            for p in new
                remotecall_wait(Core.eval,p,Main,quote
                    using LinearAlgebra,SHA
                    LinearAlgebra.BLAS.set_num_threads(1)
                    @assert VERSION==$(VERSION) "Julia versions differ"
                    for (f,h) in zip(["CtsemJulia.jl","CtsemCompute.jl"],$hashes)
                        @assert bytes2hex(sha256(read(joinpath($path,f))))==h "Backend code differs"
                        include(joinpath($path,f))
                    end
                end)
            end
        end
        isempty(ids)&&error("No compute workers")
        ids
    catch
        isempty(ids)||rmprocs(ids...)
        rethrow()
    end
end
function distribute!(ids,sd)
    n=Int(CtsemJulia.field(sd,:nsubjects))
    length(ids)<=n||error("More workers than independent subjects")
    for (j,p) in enumerate(ids)
        block=fld((j-1)*n,length(ids))+1:fld(j*n,length(ids))
        remotecall_wait(setup,p,sd,block)
    end
end
function jointfg(x,ids,chunksize)
    pending=[remotecall(evaluate,p,x,chunksize) for p in ids]
    parts=fetch.(pending)
    (value=sum(p.value for p in parts),gradient=reduce(+,p.gradient for p in parts))
end
function atomic_save(file,x)
    isempty(file)&&return
    mkpath(dirname(file));tmp=file*".tmp-"*string(getpid())
    serialize(tmp,x);mv(tmp,file;force=true)
end
function fingerprint(sd,x,options)
    io=IOBuffer();serialize(io,(sd,x,options,read(@__FILE__),read(joinpath(@__DIR__,"CtsemJulia.jl"))))
    bytes2hex(sha256(take!(io)))
end
function summary_result(x,start,fg,iterations,ok,scale=ones(length(x)))
    r=fg(x);g=r.gradient;g0=fg(start).gradient
    (estimates=x,loglik=-r.value,iterations=iterations,stage_iterations=[iterations],
     stage_objectives=[r.value],stage_maximum_absolute_gradients=[maximum(abs.(g))],
     optimizer_converged=ok,termination_code=ok ? "Converged" : "NotConverged",stopped_by="Compute",
     converged=ok,initial_maximum_absolute_gradient=maximum(abs.(g0)),
     maximum_absolute_gradient=maximum(abs.(g)),gradient_L2_norm=norm(g),
     maximum_parameter_change=maximum(abs.(x-start)),parameter_change_L2=norm(x-start),
     scale_minimum=minimum(scale),scale_maximum=maximum(scale))
end
function joint(start,sd,ids;maxiters=5000,chunksize=7,gradient_tolerance=1e-3,
    initial_raw_step=.01,show_trace=true,checkpoint="",signature="",function_relative_tolerance=0.,parameter_relative_tolerance=0.)
    distribute!(ids,sd)
    fg=x->jointfg(x,ids,chunksize)
    r=fg(start)
    expected=CtsemJulia.ctsem_negloglik(start,sd)
    expectedg=CtsemJulia.ctsem_gradient_by_subject(start,sd;threaded=false,chunksize=chunksize)
    isapprox(r.value,expected;rtol=1e-10,atol=1e-7)||error("Distributed likelihood validation failed")
    maximum(abs.(r.gradient-expectedg)./(1 .+abs.(expectedg)))<1e-7||error("Distributed gradient validation failed")
    origin=copy(start)
    scale=initial_raw_step ./ max.(abs.(r.gradient),1e-3)
    scale=clamp.(scale,1e-6,1e2)
    best=Ref(copy(start));bestf=Ref(r.value);cachex=fill(NaN,length(start));cached=Ref{Any}(nothing)
    function evaluate_u(u)
        x=origin+scale.*u
        if cached[]===nothing || !isequal(x,cachex)
            cached[]=fg(x);copyto!(cachex,x)
            if isfinite(cached[].value)&&cached[].value<bestf[]
                best[]=copy(x);bestf[]=cached[].value
                atomic_save(checkpoint,(signature=signature,estimates=best[],completed=false))
            end
        end
        cached[]
    end
    f=u->evaluate_u(u).value
    g! = (g,u)->(g .= scale .* evaluate_u(u).gradient)
    result=Optim.optimize(f,g!,zeros(length(start)),Optim.LBFGS(),Optim.Options(
        iterations=maxiters,g_abstol=max(eps(),minimum(scale)*gradient_tolerance),
        f_reltol=function_relative_tolerance,x_reltol=parameter_relative_tolerance,show_trace=show_trace))
    final=fg(best[])
    summary_result(best[],start,fg,Optim.iterations(result),maximum(abs.(final.gradient))<=gradient_tolerance,scale)
end
function one_start(x,sd,options)
    o=CtsemJulia.ctsem_optimize_auto(x,sd;options...)
    fg=z->(value=CtsemJulia.ctsem_negloglik(z,sd),gradient=CtsemJulia.ctsem_gradient_by_subject(z,sd;threaded=false))
    r=summary_result(o.estimates,x,fg,o.total_iterations,
        maximum(abs.(fg(o.estimates).gradient))<=options.gradient_tolerance,o.scale)
    merge(r,(host=gethostname(),))
end
function retryable(err)
    err isa ProcessExitedException || err isa EOFError || err isa IOError ||
    (err isa CompositeException && any(retryable,err.exceptions)) ||
    (err isa TaskFailedException && retryable(err.task.result)) ||
    (err isa RemoteException && retryable(err.captured.ex))
end
function optimize(start_r,sd_r,config;tasks="likelihood",starts=1,jitter=.1,seed=1,
    maxiters=5000,chunksize=7,gradient_tolerance=1e-3,initial_raw_step=.01,show_trace=true,
    checkpoint="",retries=0,function_relative_tolerance=0.,parameter_relative_tolerance=0.)
    start=Float64.(CtsemJulia.fvec(start_r));sd=Dict(Symbol(k)=>(string(parentmodule(typeof(v)))=="DataFrames" ? Matrix(v) : v) for (k,v) in pairs(sd_r))
    # RCall data is normalized by the R adapter before entering this function.
    CtsemJulia.check_supported(sd)
    Int(get(sd,:nstudyre,0))==0||error("Distributed ML does not support integrated study random effects")
    Int(CtsemJulia.field(sd,:priors))==0||error("Distributed ML does not support priors")
    Int(CtsemJulia.field(sd,:nindvarying))==0||error("Distributed ML does not support integrated random effects")
    options=(maxiters=Int(maxiters),chunksize=Int(chunksize),gradient_tolerance=Float64(gradient_tolerance),
        initial_raw_step=Float64(initial_raw_step),show_trace=false,threaded=false,
        function_relative_tolerance=Float64(function_relative_tolerance),parameter_relative_tolerance=Float64(parameter_relative_tolerance))
    signature=fingerprint(sd,start,(options,tasks,starts,jitter,seed))
    if !isempty(checkpoint)&&isfile(checkpoint)
        saved=deserialize(checkpoint)
        saved.signature==signature||error("Checkpoint belongs to different data/model/settings")
        saved.completed&&return saved.result
        tasks=="likelihood"&&(start=copy(saved.estimates))
    end
    for attempt in 0:retries
        ids=Int[]
        try
            ids=connect(config)
            result=if tasks=="likelihood"
                joint(start,sd,ids;maxiters=maxiters,chunksize=chunksize,gradient_tolerance=gradient_tolerance,
                    initial_raw_step=initial_raw_step,show_trace=show_trace,checkpoint=checkpoint,signature=signature,
                    function_relative_tolerance=function_relative_tolerance,parameter_relative_tolerance=parameter_relative_tolerance)
            elseif tasks=="starts"
                vectors=[start .+ (i==1 ? zeros(length(start)) : jitter.*max.(abs.(start),.1).*randn(MersenneTwister(seed+i),length(start))) for i in 1:starts]
                values=Vector{Any}(undef,length(vectors))
                pool=Channel{Int}(length(ids));foreach(p->put!(pool,p),ids)
                @sync for i in eachindex(vectors)
                    @async begin
                        file=isempty(checkpoint) ? "" : checkpoint*".start-"*string(i)
                        cached=(!isempty(file)&&isfile(file)) ? deserialize(file) : nothing
                        if cached!==nothing
                            cached.signature==signature||error("Start checkpoint mismatch")
                            values[i]=cached.result
                        else
                            p=take!(pool)
                            try
                                values[i]=remotecall_fetch(one_start,p,vectors[i],sd,options)
                                atomic_save(file,(signature=signature,result=values[i]))
                            finally
                                put!(pool,p)
                            end
                        end
                    end
                end
                valid=findall(r->r.converged&&isfinite(r.loglik),values)
                isempty(valid)&&error("No distributed start converged; no successful fit selected")
                values[valid[argmax([values[i].loglik for i in valid])]]
            else
                error("Unsupported compute task: "*tasks)
            end
            result=merge(result,(compute=(hosts=[remotecall_fetch(gethostname,p) for p in ids],tasks=tasks,attempt=attempt),))
            atomic_save(checkpoint,(signature=signature,completed=true,result=result))
            return result
        catch err
            # Only network/worker loss is retryable; model and validation errors are not.
            (attempt<retries&&retryable(err))||rethrow()
            println("Compute connection lost: restarting workers and current task (attempt ",attempt+1,")")
            if !isempty(checkpoint)&&isfile(checkpoint)
                saved=deserialize(checkpoint);!saved.completed&&tasks=="likelihood"&&(start=copy(saved.estimates))
            end
        finally
            isempty(ids)||try rmprocs(ids...) catch end
        end
    end
end
end
