module CtsemJulia

using LinearAlgebra
using ForwardDiff
using Optim
using Random
using Printf
import AdvancedHMC
import Distributed

export ctsem_loglik, ctsem_negloglik, ctsem_gradient, ctsem_optimize,
       ctsem_loglik_by_subject, ctsem_gradient_by_subject,
       ctsem_re_loglik_by_subject, ctsem_re_logposterior_and_gradient,
       ctsem_hessian_by_subject, ctsem_vcov_by_subject,
       ctsem_optimize_by_subject, ctsem_optimize_auto, ctsem_optimize_staged,
       ctsem_optimize_scaled, scaled_minimizer,
       ctsem_optimize_whitened, whitened_minimizer,
       ctsem_logprior, ctsem_logposterior_and_gradient, ctsem_nuts,
       ctsem_nuts_parallel, ctsem_optimize_multistart,
       check_supported, tform

const LOG2PI = log(2pi)

"Create a compact Optim callback that updates one terminal line."
function optimizer_progress_callback(enabled::Bool)
    started = time()
    previous_width = Ref(0)
    function callback(state)
        enabled || return false
        iteration = hasproperty(state, :iteration) ? state.iteration :
            state.pseudo_iteration
        objective = hasproperty(state, :value) ? state.value : state.f_x
        gradient = hasproperty(state, :g_norm) ? state.g_norm :
            maximum(abs, state.g_x)
        line = @sprintf(
            "Iter %6d | objective %.8e | gradient %.6e | time %.1f s",
            iteration, objective, gradient, time() - started)
        padding = max(0, previous_width[] - ncodeunits(line))
        print('\r', line, repeat(" ", padding))
        flush(stdout)
        previous_width[] = ncodeunits(line)
        false
    end
    callback
end

"Validate and unpack the first supported CoTiMA study-RE parameterization."
function prepare_study_re(x, sd, d)
    q = Int(field(sd, :nstudyre))
    K = Int(field(sd, :nstudies))
    q > 0 || error("The study-RE model contains no random effects")
    types = ivec(field(sd, :studyre_type))
    all(==(1), types) || error(
        "Julia study-RE fitting currently supports DRIFT/raw type-1 effects only")
    indices = ivec(field(sd, :studyre_index))
    all(i -> 1 <= i <= d.nparams, indices) ||
        error("A study-RE raw parameter index is out of range")
    study_id = ivec(field(sd, :study_id))
    length(study_id) == d.nsubjects || error(
        "study_id must contain one entry per subject")
    all(s -> 1 <= s <= K, study_id) || error("study_id is out of range")
    B = d.nparams + d.ntipredeffects
    expected = B + q + q * K
    length(x) == expected || error(
        "Study-RE parameter length $(length(x)); expected $expected " *
        "($B ctsem + $q log-SD + $(q*K) standardized effects)")
    base = @view x[1:B]
    logsd = @view x[B+1:B+q]
    z = reshape(@view(x[B+q+1:end]), q, K)
    tau = exp.(logsd)
    effects = tau .* z
    (
        q=q, K=K, B=B, base=base, logsd=logsd, z=z, tau=tau,
        effects=effects, indices=indices, study_id=study_id,
        sd_prior=fvec(field(sd, :studyre_sd_prior)),
        use_hyperpriors=Int(field(sd, :studyre_use_hyperpriors)) == 1
    )
end

"Construct the ordinary population and TI-effect contribution for a subject."
function subject_raw_from_base(base, d, si)
    rawpop = @view base[1:d.nparams]
    tipar = @view base[d.nparams+1:end]
    localraw = copy(rawpop)
    for j in 1:d.ntipred, p in 1:d.nparams
        ix = d.tie_setup[p, j]
        ix > 0 && (localraw[p] += tipar[ix] * d.tipreds[si, j])
    end
    localraw
end

"Likelihood for independent noncentered study effects on raw DRIFT parameters."
function ctsem_re_loglik_by_subject(rawest_r, sd; threaded::Bool=false)
    check_supported(sd)
    x = fvec(rawest_r)
    d = prepare_subject_data(sd)
    re = prepare_study_re(x, sd, d)
    function one_subject(si)
        localraw = subject_raw_from_base(re.base, d, si)
        s = re.study_id[si]
        for r in 1:re.q
            localraw[re.indices[r]] += re.effects[r, s]
        end
        subject_loglik(localraw, d, d.firstrow[si], d.lastrow[si])
    end
    if threaded && Threads.nthreads() > 1
        thread_ll = zeros(eltype(x), Threads.maxthreadid())
        Threads.@threads :static for si in 1:d.nsubjects
            thread_ll[Threads.threadid()] += one_subject(si)
        end
        sum(thread_ll)
    else
        sum(one_subject(si) for si in 1:d.nsubjects)
    end
end

"Log posterior and gradient for the first noncentered independent study-RE slice."
function ctsem_re_logposterior_and_gradient(rawest_r, sd;
                                            chunksize::Integer=7,
                                            threaded::Bool=false)
    check_supported(sd)
    x = fvec(rawest_r)
    d = prepare_subject_data(sd)
    re = prepare_study_re(x, sd, d)
    re.use_hyperpriors || error(
        "The initial Julia study-RE implementation supports Bayesian hyperpriors only")
    length(re.sd_prior) == re.q || error("studyre_sd_prior has the wrong length")
    all(>(0), re.sd_prior) || error("studyre_sd_prior must be positive")

    function add_subject!(gbase, geffect, si)
        localraw = subject_raw_from_base(re.base, d, si)
        s = re.study_id[si]
        for r in 1:re.q
            localraw[re.indices[r]] += re.effects[r, s]
        end
        f = y -> -subject_loglik(y, d, d.firstrow[si], d.lastrow[si])
        cs = min(chunksize, length(localraw))
        cfg = ForwardDiff.GradientConfig(f, localraw, ForwardDiff.Chunk{cs}())
        gs = ForwardDiff.gradient(f, localraw, cfg)
        gbase[1:d.nparams] .+= gs
        for j in 1:d.ntipred, p in 1:d.nparams
            ix = d.tie_setup[p, j]
            ix > 0 && (gbase[d.nparams + ix] += gs[p] * d.tipreds[si, j])
        end
        for r in 1:re.q
            geffect[r, s] += gs[re.indices[r]]
        end
        nothing
    end

    if threaded && Threads.nthreads() > 1
        gb = [zeros(eltype(x), re.B) for _ in 1:Threads.maxthreadid()]
        ge = [zeros(eltype(x), re.q, re.K) for _ in 1:Threads.maxthreadid()]
        Threads.@threads :static for si in 1:d.nsubjects
            tid = Threads.threadid()
            add_subject!(gb[tid], ge[tid], si)
        end
        nll_base = reduce(+, gb)
        nll_effect = reduce(+, ge)
    else
        nll_base = zeros(eltype(x), re.B)
        nll_effect = zeros(eltype(x), re.q, re.K)
        for si in 1:d.nsubjects
            add_subject!(nll_base, nll_effect, si)
        end
    end

    ll = ctsem_re_loglik_by_subject(x, sd; threaded=threaded)
    lp = ctsem_logprior(re.base, sd)
    lp -= (sum(abs2, re.z) + length(re.z) * LOG2PI) / 2
    for r in 1:re.q
        sigma = re.sd_prior[r]
        lp += -(re.tau[r] / sigma)^2 / 2 - log(sigma) - LOG2PI / 2 + re.logsd[r]
    end

    g = zeros(eltype(x), length(x))
    g[1:re.B] .= -nll_base
    nparams = d.nparams
    priormod = field(sd, :priormod) / Int(field(sd, :nsubsets))
    scale = field(sd, :tipredeffectscale)
    g[1:nparams] .-= priormod .* @view(re.base[1:nparams])
    d.ntipredeffects > 0 && (g[nparams+1:re.B] .-=
        (priormod / scale^2) .* @view(re.base[nparams+1:end]))
    loglik_effect = -nll_effect
    for r in 1:re.q
        g[re.B + r] = sum(@view(loglik_effect[r, :]) .* @view(re.z[r, :])) *
            re.tau[r] + 1 - (re.tau[r] / re.sd_prior[r])^2
    end
    zoffset = re.B + re.q
    for s in 1:re.K, r in 1:re.q
        pos = zoffset + (s - 1) * re.q + r
        g[pos] = loglik_effect[r, s] * re.tau[r] - re.z[r, s]
    end
    (ll + lp, g)
end

"Access a field from a Julia Dict, NamedTuple, or struct converted by JuliaCall."
function field(x, name::Symbol)
    if x isa AbstractDict
        haskey(x, name) && return x[name]
        haskey(x, String(name)) && return x[String(name)]
    end
    hasproperty(x, name) && return getproperty(x, name)
    error("standata has no field $(name)")
end

ivec(x::Number) = [Int(x)]
ivec(x) = Int.(vec(x))

fvec(x::Number) = [x]
fvec(x) = collect(vec(x))

# Avoid `abs`/`max` here: although value-stable, their branch derivatives at
# exactly zero give ForwardDiff the wrong sign. This form is equally stable and
# has the correct softplus derivative of 1/2 at zero.
softplus(x) = x > zero(x) ? x + log1p(exp(-x)) : log1p(exp(x))
logistic(x) = inv(one(x) + exp(-x))

function tform(parin, transform::Integer, scale, meanscale, shift, innershift)
    p = parin
    meanscale != 1 && (p *= meanscale)
    innershift != 0 && (p += innershift)
    transform == 1  && (p = softplus(p))
    transform == 2  && (p = exp(p))
    transform == 3  && (p = logistic(p))
    transform == 4  && (p = p^3)
    transform == 5  && (p = log1p(p))
    transform == 50 && (p = meanscale)
    transform == 51 && (p = logistic(p))
    transform == 52 && (p = exp(p))
    transform == 53 && (p = logistic(p) - exp(2p) / (1 + exp(p))^2)
    transform == 54 && (p = 3p^2)
    transform == 55 && (p = inv(1 + p))
    scale != 1 && (p *= scale)
    transform < 49 && shift != 0 && (p += shift)
    p
end

"Convert a ctsem SD/correlation-square-root parameter matrix to covariance."
function constraincorsqrt1(mat::AbstractMatrix)
    d = size(mat, 1)
    T = eltype(mat)
    out = zeros(T, d, d)
    for i in 1:d
        s = T(1e-5)
        ss = T(1e-5)
        for j in 1:d
            if j > i
                ss += mat[j, i]^2
                s += mat[j, i]
            elseif j < i
                ss += mat[i, j]^2
                s += mat[i, j]
            end
        end
        r1 = sqrt(ss)
        r3 = abs(s) / r1 - 1
        r4 = sqrt(softplus(2 * (abs(s) - s - 1) - 4))
        r = sqrt(ss + (r4 * r3 + 1) * r4 + 1)
        for j in 1:d
            j > i && (out[i, j] = mat[j, i] / r)
            j < i && (out[i, j] = mat[i, j] / r)
        end
        out[i, i] = sqrt(1 - sum(abs2, @view(out[i, :])) + T(1e-5))
    end
    out
end

function sdcovsqrt2cov(mat::AbstractMatrix, choleskymats::Integer)
    isempty(mat) && return copy(mat)
    choleskymats < 1 ? begin
        c = constraincorsqrt1(mat)
        l = Diagonal(diag(mat)) * c
        l * l'
    end : mat * mat'
end

function makesym(a::AbstractMatrix)
    out = copy(a)
    for j in axes(out, 2), i in axes(out, 1)
        if i == j
            out[i, j] = a[i, j] + 1e-10
        elseif i < j
            out[i, j] = a[j, i]
        end
    end
    out
end

"Solve A*P + P*A' + Q = 0; equivalent to ctsem's ksolve for symmetric P."
function ksolve(A::AbstractMatrix, Q::AbstractMatrix)
    d = size(A, 1)
    T = promote_type(eltype(A), eltype(Q))
    # Vectorised Lyapunov equation. Column-major vec convention:
    # vec(A*P + P*A') = (I kron A + A kron I) vec(P).
    K = kron(Matrix{T}(I, d, d), A) + kron(A, Matrix{T}(I, d, d))
    p = -(K \ vec(Q))
    Matrix(Symmetric(reshape(p, d, d)))
end

function expm_subsets(m::AbstractMatrix, subsets::AbstractMatrix{<:Integer})
    n = size(m, 1)
    out = zeros(eltype(m), n, n)
    for r in axes(subsets, 1)
        ix = [j for j in @view(subsets[r, :]) if j != 0]
        isempty(ix) && continue
        if length(ix) == 1
            out[ix[1], ix[1]] = exp(m[ix[1], ix[1]])
        else
            out[ix, ix] = matrix_exp(m[ix, ix])
        end
    end
    out
end

# Base/BLAS matrix exponential for the already validated Float64 objective.
matrix_exp(a::Matrix{Float64}) = exp(a)

"""Generic [13/13] Padé matrix exponential for automatic-differentiation types."""
function matrix_exp(a::AbstractMatrix)
    # Fixed scaling avoids differentiating a discrete norm-based branch. Eight
    # squarings are conservative for the small ctsem drift blocks in this fit.
    s = 8
    A = a / (2^s)
    n = size(A, 1)
    T = eltype(A)
    Id = Matrix{T}(I, n, n)
    A2 = A * A
    A4 = A2 * A2
    A6 = A4 * A2
    b = (64764752532480000.0, 32382376266240000.0,
         7771770303897600.0, 1187353796428800.0,
         129060195264000.0, 10559470521600.0,
         670442572800.0, 33522128640.0, 1323241920.0,
         40840800.0, 960960.0, 16380.0, 182.0, 1.0)
    U = A * (A6 * (b[14] * A6 + b[12] * A4 + b[10] * A2) +
             b[8] * A6 + b[6] * A4 + b[4] * A2 + b[2] * Id)
    V = A6 * (b[13] * A6 + b[11] * A4 + b[9] * A2) +
        b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * Id
    R = (V - U) \ (V + U)
    for _ in 1:s
        R = R * R
    end
    R
end

function parvectform(which::AbstractVector{<:Integer}, rawpar, when::Integer,
                     ms, mval, subi::Integer)
    out = similar(rawpar, length(which))
    for (oi, pindex) in pairs(which)
        found = false
        for r in axes(ms, 1)
            if !found && (ms[r, 8] == when || ms[r, 8] == 100) && ms[r, 3] == pindex
                if subi == 0 || (ms[r, 3] > 0 &&
                    (ms[r, 5] > 0 || ms[r, 6] > 0 || ms[r, 8] > 0))
                    out[oi] = tform(rawpar[pindex], ms[r, 4], mval[r, 2],
                                    mval[r, 3], mval[r, 4], mval[r, 6])
                end
                found = true
            end
        end
        found || error("No matsetup row found for parameter $pindex")
    end
    out
end

function mcalc!(out, tfpars, states, when, matrix_number, ms, mval, subi)
    changed = false
    for r in axes(ms, 1)
        if ms[r, 7] == matrix_number &&
           (subi == 0 || (ms[r, 3] > 0 &&
             (ms[r, 5] > 0 || ms[r, 6] > 0 || ms[r, 8] > 0)))
            whenyes = ms[r, 8] == 100 || any(==(ms[r, 8]), when)
            if whenyes
                changed = true
                rr, cc, pi = ms[r, 1], ms[r, 2], ms[r, 3]
                if pi > 0 && ms[r, 8] == 0
                    out[rr, cc] = tfpars[pi]
                elseif pi > 0
                    out[rr, cc] = tform(states[pi], ms[r, 4], mval[r, 2],
                                        mval[r, 3], mval[r, 4], mval[r, 6])
                else
                    out[rr, cc] = mval[r, 1]
                end
            end
        end
    end
    changed
end

function blank_matrix(matrixdims, i, T=Float64)
    fill(T(NaN), Int(matrixdims[i, 1]), Int(matrixdims[i, 2]))
end

function update_parameter_vector!(indparams, rawparams, row, whenvecp, ms, mval, subi)
    which = findall(!=(0), @view whenvecp[row, :])
    vals = parvectform(which, rawparams, 0, ms, mval, subi)
    indparams[which] .= vals
end

function check_supported(sd; throw_error=true)
    problems = String[]
    expected = Dict(
        :continuoustime => 1, :nindvarying => 0, :intoverstates => 1,
        :savescores => 0, :savesubjectmatrices => 0,
        :dokalman => 1, :nsubsets => 1, :ntdpred => 0)
    for (k, v) in expected
        Int(field(sd, k)) == v || push!(problems, "$k must equal $v")
    end
    Int(field(sd, :priors)) in (0, 1) ||
        push!(problems, "priors must equal 0 or 1")
    all(==(0), ivec(field(sd, :manifesttype))) ||
        push!(problems, "only continuous manifest variables are supported")
    Int(field(sd, :nmissingtipreds)) == 0 ||
        push!(problems, "missing time-independent predictors are not supported")
    all(==(0), ivec(field(sd, :statedep))) ||
        push!(problems, "state-dependent parameters are not yet supported")
    Int(field(sd, :nJAxfinite)) == 0 ||
        push!(problems, "finite-difference drift Jacobians are not yet supported")
    Int(field(sd, :nJyfinite)) == 0 ||
        push!(problems, "finite-difference measurement Jacobians are not yet supported")
    if throw_error && !isempty(problems)
        error("Unsupported ctsem branch(es):\n  - " * join(problems, "\n  - "))
    end
    problems
end

function mvnormal_logpdf_chol(y, mu, cov)
    C = cholesky(Symmetric(cov); check=true)
    z = C.L \ (y - mu)
    -0.5 * (length(y) * LOG2PI + 2sum(log, diag(C.L)) + dot(z, z))
end

"""
    ctsem_loglik(rawest, sd; rows_limit=0, return_rows=false)

Specialised Julia port of the linear-Gaussian ctsem likelihood represented by
the supplied CoTiMA RI fit. `rawest` must contain `nparams` population means
followed by `ntipredeffects` TI-predictor effects, exactly as in ctsem's
optimisation vector. Set `rows_limit` to a small number for initial debugging;
zero uses every row.
"""
function ctsem_loglik(rawest_r, sd; rows_limit::Integer=0, return_rows::Bool=false)
    check_supported(sd)
    rawest = fvec(rawest_r)
    nparams = Int(field(sd, :nparams))
    ntie = Int(field(sd, :ntipredeffects))
    length(rawest) == nparams + ntie ||
        error("rawest length $(length(rawest)); expected $(nparams + ntie)")

    rawpopmeans = rawest[1:nparams]
    tipar = rawest[nparams+1:end]
    ms = Int.(field(sd, :matsetup))
    mval = field(sd, :matvalues)
    matrixdims = Int.(field(sd, :matrixdims))
    whenmat = Int.(field(sd, :whenmat))
    whenvecp = Int.(field(sd, :whenvecp))
    subsets = Int.(field(sd, :JAxsubsets))
    tie_setup = Int.(field(sd, :TIPREDEFFECTsetup))
    tipreds = field(sd, :tipredsdata)
    Y = field(sd, :Y)
    times = fvec(field(sd, :time))
    subjects = ivec(field(sd, :subject))
    dokalmanrows = ivec(field(sd, :dokalmanrows))
    nobs = ivec(field(sd, :nobs_y))
    whichobs = Int.(field(sd, :whichobs_y))
    ncont = ivec(field(sd, :ncont_y))
    whichcont = Int.(field(sd, :whichcont_y))
    derrind = ivec(field(sd, :derrind))
    maxtimestep = Float64(field(sd, :maxtimestep))
    choleskymats = Int(field(sd, :choleskymats))
    nlatent = Int(field(sd, :nlatent))
    nlatentpop = Int(field(sd, :nlatentpop))
    nd = length(subjects)
    lastrow = rows_limit > 0 ? min(rows_limit, nd) : nd

    # Design matrix from the compact ctsem indices.
    TIP = zeros(eltype(rawest), nparams, Int(field(sd, :ntipred)))
    for j in axes(TIP, 2), i in axes(TIP, 1)
        ix = tie_setup[i, j]
        ix > 0 && (TIP[i, j] = tipar[ix])
    end

    mats = Dict{Int, Matrix{eltype(rawest)}}()
    for i in (1,2,3,4,5,6,7,8,9,10,31,32,33,21,22,52,53,54)
        mats[i] = blank_matrix(matrixdims, i, eltype(rawest))
    end

    # Population initialisation (Stan rowx == 0).
    indparams = similar(rawpopmeans)
    update_parameter_vector!(indparams, rawpopmeans, 1, whenvecp, ms, mval, 0)
    dummy_state = fill(eltype(rawest)(-999), nlatentpop)
    for (mi, whens) in ((1,[0]), (10,[0,1]), (1,[0,1]), (8,[0,1]),
                        (2,[0]), (3,[0]), (4,[0]), (5,[0]), (6,[0]),
                        (7,[0]), (8,[0]), (9,[0]), (52,[0]),
                        (53,[0]), (54,[0]))
        mcalc!(mats[mi], indparams, dummy_state, whens, mi, ms, mval, 0)
        mi == 1 && (dummy_state .= vec(mats[1]))
    end

    llrows = zeros(eltype(rawest), lastrow)
    prevrow = 0
    prevdt = NaN
    current_subject = 0
    state = zeros(eltype(rawest), nlatentpop)
    etacov = zeros(eltype(rawest), nlatentpop, nlatentpop)
    diffusioncov = zeros(eltype(rawest), nlatent, nlatent)
    asymdiff = zeros(eltype(rawest), nlatent, nlatent)

    for row in 1:lastrow
        dokalmanrows[row] == 0 && continue
        si = subjects[row]
        new_subject = si != current_subject

        if new_subject
            current_subject = si
            rawparams = rawpopmeans + TIP * vec(@view tipreds[si, :])
            update_parameter_vector!(indparams, rawparams, 2, whenvecp, ms, mval, si)
            # Reproduce ctsem's subject-initialisation order. This matters when
            # one matrix refers to a parameter represented in PARS/T0MEANS.
            whenmat[1, 5] >= 1 && mcalc!(mats[1], indparams, state, [0], 1, ms, mval, si)
            state .= vec(mats[1])
            for (mi, whens) in ((10,[0,1]), (1,[0,1]), (8,[0,1]),
                                (2,[0]), (3,[0]), (4,[0]), (5,[0]),
                                (6,[0]), (7,[0]), (8,[0]), (9,[0]),
                                (52,[0]), (53,[0]), (54,[0]))
                whenmat[mi, 5] > 0 &&
                    mcalc!(mats[mi], indparams, state, whens, mi, ms, mval, si)
            end
            state .= vec(mats[1])
            t0cov = sdcovsqrt2cov(mats[8], choleskymats)
            etacov .= t0cov
        else
            dt = times[row] - times[prevrow]
            dtsmall = dt / ceil(dt / maxtimestep)
            elapsed = zero(dt)
            while elapsed < dt - 1e-10
                elapsed += dtsmall
                drift = mats[3]
                cint = vec(mats[7])
                JAx = mats[52]
                diffusioncov[derrind, derrind] =
                    sdcovsqrt2cov(mats[4][derrind, derrind], choleskymats)
                eJAx = expm_subsets(JAx * dtsmall, subsets)
                if isnan(prevdt) || dt != prevdt || elapsed == dtsmall
                    asymdiff[derrind, derrind] =
                        ksolve(JAx[derrind, derrind], diffusioncov[derrind, derrind])
                end
                discrete_diff = asymdiff[derrind, derrind] -
                    eJAx[derrind, derrind] * asymdiff[derrind, derrind] *
                    eJAx[derrind, derrind]'
                f0 = drift * state + cint
                affine = f0 - JAx * state
                dcint = JAx \ ((eJAx - I) * affine)
                state .= eJAx * state + dcint
                etacov .= eJAx * makesym(etacov) * eJAx'
                etacov[derrind, derrind] .+= discrete_diff
            end
            prevdt = dt
        end

        od = nobs[row] == 0 ? Int[] : vec(whichobs[row, 1:nobs[row]])
        o0 = ncont[row] == 0 ? Int[] : vec(whichcont[row, 1:ncont[row]])
        if !isempty(od)
            lambda = mats[2]
            manifestmeans = vec(mats[6])
            manifestcov = sdcovsqrt2cov(mats[5], choleskymats)
            Jy = mats[54]
            sy = lambda * state + manifestmeans
            ycov = Jy * makesym(etacov) * Jy' + manifestcov
            err = Y[row, od] - sy[od]
            K = (etacov * Jy[od, :]') / Symmetric(ycov[od, od])
            etacov .-= K * Jy[od, :] * etacov
            state .+= K * err
            !isempty(o0) && (llrows[row] =
                mvnormal_logpdf_chol(Y[row, o0], sy[o0], ycov[o0, o0]))
        end
        prevrow = row
    end

    result = (loglik=sum(llrows), rows_used=lastrow,
              nonzero_rows=count(!iszero, llrows))
    return_rows ? merge(result, (llrow=llrows,)) : result
end

ctsem_negloglik(rawest, sd) = -ctsem_loglik(rawest, sd).loglik

"Prepare immutable data views and subject row ranges for score decomposition."
function prepare_subject_data(sd)
    subjects = ivec(field(sd, :subject))
    nsubjects = Int(field(sd, :nsubjects))
    firstrow = zeros(Int, nsubjects)
    lastrow = zeros(Int, nsubjects)
    for r in eachindex(subjects)
        si = subjects[r]
        firstrow[si] == 0 && (firstrow[si] = r)
        lastrow[si] = r
    end
    any(==(0), firstrow) && error("At least one subject has no active data row")
    ntipred = Int(field(sd, :ntipred))
    tipreds = ntipred == 0 ? zeros(nsubjects, 0) : field(sd, :tipredsdata)
    subsets_raw = field(sd, :JAxsubsets)
    subsets = subsets_raw isa Number ? reshape([Int(subsets_raw)], 1, 1) :
        Int.(subsets_raw)
    (
        nparams=Int(field(sd, :nparams)),
        ntipred=Int(field(sd, :ntipred)),
        ntipredeffects=Int(field(sd, :ntipredeffects)),
        nlatent=Int(field(sd, :nlatent)),
        nlatentpop=Int(field(sd, :nlatentpop)),
        nsubjects=nsubjects,
        choleskymats=Int(field(sd, :choleskymats)),
        maxtimestep=field(sd, :maxtimestep),
        ms=Int.(field(sd, :matsetup)),
        mval=field(sd, :matvalues),
        matrixdims=Int.(field(sd, :matrixdims)),
        whenvecp=Int.(field(sd, :whenvecp)),
        subsets=subsets,
        tie_setup=Int.(field(sd, :TIPREDEFFECTsetup)),
        tipreds=tipreds,
        Y=field(sd, :Y),
        times=fvec(field(sd, :time)),
        nobs=ivec(field(sd, :nobs_y)),
        whichobs=Int.(field(sd, :whichobs_y)),
        ncont=ivec(field(sd, :ncont_y)),
        whichcont=Int.(field(sd, :whichcont_y)),
        derrind=ivec(field(sd, :derrind)),
        firstrow=firstrow,
        lastrow=lastrow
    )
end

"Construct all matrices directly from one subject's 21 raw parameters."
function subject_matrices(rawparams, d)
    T = eltype(rawparams)
    mats = Dict{Int, Matrix{T}}()
    for i in (1,2,3,4,5,6,7,8,9,10,31,32,33,21,22,52,53,54)
        mats[i] = blank_matrix(d.matrixdims, i, T)
    end
    indparams = similar(rawparams)
    update_parameter_vector!(indparams, rawparams, 1,
                             d.whenvecp, d.ms, d.mval, 0)
    state = fill(T(-999), d.nlatentpop)
    for (mi, whens) in ((1,[0]), (10,[0,1]), (1,[0,1]), (8,[0,1]),
                        (2,[0]), (3,[0]), (4,[0]), (5,[0]), (6,[0]),
                        (7,[0]), (8,[0]), (9,[0]), (52,[0]),
                        (53,[0]), (54,[0]))
        mcalc!(mats[mi], indparams, state, whens, mi, d.ms, d.mval, 0)
        mi == 1 && (state .= vec(mats[1]))
    end
    mats
end

"Likelihood contribution of one independent subject on its local 21-par scale."
function subject_loglik(rawparams, d, first::Integer, last::Integer)
    mats = subject_matrices(rawparams, d)
    T = eltype(rawparams)
    state = copy(vec(mats[1]))
    etacov = sdcovsqrt2cov(mats[8], d.choleskymats)
    drift = mats[3]
    cint = vec(mats[7])
    JAx = mats[52]
    lambda = mats[2]
    manifestmeans = vec(mats[6])
    manifestcov = sdcovsqrt2cov(mats[5], d.choleskymats)
    Jy = mats[54]
    diffusioncov = zeros(T, d.nlatent, d.nlatent)
    diffusioncov[d.derrind, d.derrind] =
        sdcovsqrt2cov(mats[4][d.derrind, d.derrind], d.choleskymats)
    asymdiff = zeros(T, d.nlatent, d.nlatent)
    asymdiff[d.derrind, d.derrind] =
        ksolve(JAx[d.derrind, d.derrind],
               diffusioncov[d.derrind, d.derrind])
    ll = zero(T)

    for row in first:last
        if row > first
            dt = d.times[row] - d.times[row - 1]
            dtsmall = dt / ceil(dt / d.maxtimestep)
            elapsed = zero(dt)
            while elapsed < dt - 1e-10
                elapsed += dtsmall
                eJAx = expm_subsets(JAx * dtsmall, d.subsets)
                discrete_diff = asymdiff[d.derrind, d.derrind] -
                    eJAx[d.derrind, d.derrind] *
                    asymdiff[d.derrind, d.derrind] *
                    eJAx[d.derrind, d.derrind]'
                f0 = drift * state + cint
                affine = f0 - JAx * state
                dcint = JAx \ ((eJAx - I) * affine)
                state .= eJAx * state + dcint
                etacov .= eJAx * makesym(etacov) * eJAx'
                etacov[d.derrind, d.derrind] .+= discrete_diff
            end
        end

        od = d.nobs[row] == 0 ? Int[] :
             vec(d.whichobs[row, 1:d.nobs[row]])
        o0 = d.ncont[row] == 0 ? Int[] :
             vec(d.whichcont[row, 1:d.ncont[row]])
        if !isempty(od)
            sy = lambda * state + manifestmeans
            ycov = Jy * makesym(etacov) * Jy' + manifestcov
            err = d.Y[row, od] - sy[od]
            K = (etacov * Jy[od, :]') / Symmetric(ycov[od, od])
            etacov .-= K * Jy[od, :] * etacov
            state .+= K * err
            !isempty(o0) && (ll += mvnormal_logpdf_chol(
                d.Y[row, o0], sy[o0], ycov[o0, o0]))
        end
    end
    ll
end

"Sum the independently evaluated subject likelihoods (validation function)."
function ctsem_loglik_by_subject(rawest_r, sd; threaded::Bool=false)
    check_supported(sd)
    x = fvec(rawest_r)
    d = prepare_subject_data(sd)
    rawpop = x[1:d.nparams]
    tipar = x[d.nparams+1:end]
    TIP = zeros(eltype(x), d.nparams, d.ntipred)
    for j in axes(TIP, 2), p in axes(TIP, 1)
        ix = d.tie_setup[p, j]
        ix > 0 && (TIP[p, j] = tipar[ix])
    end
    if threaded && Threads.nthreads() > 1
        # JuliaCall may execute on an interactive-pool thread whose global
        # thread ID is larger than `Threads.nthreads()` (the default-pool
        # count). Size buffers by the maximum possible global thread ID.
        thread_ll = zeros(eltype(x), Threads.maxthreadid())
        Threads.@threads :static for si in 1:d.nsubjects
            tid = Threads.threadid()
            localraw = rawpop + TIP * vec(@view d.tipreds[si, :])
            thread_ll[tid] += subject_loglik(
                localraw, d, d.firstrow[si], d.lastrow[si])
        end
        return sum(thread_ll)
    else
        ll = zero(eltype(x))
        for si in 1:d.nsubjects
            localraw = rawpop + TIP * vec(@view d.tipreds[si, :])
            ll += subject_loglik(
                localraw, d, d.firstrow[si], d.lastrow[si])
        end
        return ll
    end
end

"""
    ctsem_gradient_by_subject(rawest, sd; chunksize=7, progress_every=5000)

Differentiate each independent likelihood contribution on its 21-parameter
local scale and map the scores analytically to the 21 + 630 global parameters.
Returns the gradient of the negative log-likelihood.
"""
function ctsem_gradient_by_subject(rawest_r, sd; chunksize::Integer=7,
                                   progress_every::Integer=5000,
                                   threaded::Bool=false)
    check_supported(sd)
    x = fvec(rawest_r)
    d = prepare_subject_data(sd)
    rawpop = x[1:d.nparams]
    tipar = x[d.nparams+1:end]
    TIP = zeros(eltype(x), d.nparams, d.ntipred)
    for j in axes(TIP, 2), p in axes(TIP, 1)
        ix = d.tie_setup[p, j]
        ix > 0 && (TIP[p, j] = tipar[ix])
    end
    function add_subject_score!(g, si)
        z = vec(@view d.tipreds[si, :])
        localraw = rawpop + TIP * z
        f = q -> -subject_loglik(q, d, d.firstrow[si], d.lastrow[si])
        local_chunksize = min(chunksize, length(localraw))
        cfg = ForwardDiff.GradientConfig(f, localraw,
                  ForwardDiff.Chunk{local_chunksize}())
        gs = ForwardDiff.gradient(f, localraw, cfg)
        g[1:d.nparams] .+= gs
        for j in 1:d.ntipred, p in 1:d.nparams
            ix = d.tie_setup[p, j]
            ix > 0 && (g[d.nparams + ix] += gs[p] * z[j])
        end
        nothing
    end

    if threaded && Threads.nthreads() > 1
        thread_g = [zeros(eltype(x), length(x))
                    for _ in 1:Threads.maxthreadid()]
        completed = Threads.Atomic{Int}(0)
        Threads.@threads :static for si in 1:d.nsubjects
            add_subject_score!(thread_g[Threads.threadid()], si)
            done = Threads.atomic_add!(completed, 1) + 1
            progress_every > 0 && done % progress_every == 0 &&
                println("subject-score progress: ", done, "/", d.nsubjects)
        end
        g = zeros(eltype(x), length(x))
        for gt in thread_g
            g .+= gt
        end
        return g
    else
        g = zeros(eltype(x), length(x))
        for si in 1:d.nsubjects
            add_subject_score!(g, si)
            progress_every > 0 && si % progress_every == 0 &&
                println("subject-score progress: ", si, "/", d.nsubjects)
        end
        return g
    end
end

"""
    ctsem_logprior(rawest, sd)

The ordinary ctsem prior for the currently supported fixed-effects branch:
standard normal priors on raw population parameters and normal priors with
`tipredeffectscale` on raw time-independent-predictor effects. Constants are
included, although ctsem's Stan expression omits the scale Jacobian for the
predictor effects (which is constant in the sampled parameters).
"""
function ctsem_logprior(rawest_r, sd)
    x = fvec(rawest_r)
    nparams = Int(field(sd, :nparams))
    ntie = Int(field(sd, :ntipredeffects))
    length(x) == nparams + ntie || error(
        "Raw parameter length $(length(x)) does not equal nparams + " *
        "ntipredeffects ($(nparams + ntie)).")
    all(==(0), ivec(field(sd, :laplaceprior))) || error(
        "Laplace population priors are not yet supported by Julia NUTS.")
    Int(field(sd, :laplacetipreds)) == 0 || error(
        "Laplace TI-predictor priors are not yet supported by Julia NUTS.")
    Int(field(sd, :laplaceprioronly)) == 0 || error(
        "laplaceprioronly is not supported by Julia NUTS.")
    priormod = field(sd, :priormod) / Int(field(sd, :nsubsets))
    scale = field(sd, :tipredeffectscale)
    scale > 0 || error("tipredeffectscale must be positive")
    lp = -priormod * (sum(abs2, @view(x[1:nparams])) / 2 +
                      nparams * LOG2PI / 2)
    if ntie > 0
        tie = @view x[nparams+1:end]
        lp -= priormod * (sum(abs2, tie) / (2scale^2) + ntie * LOG2PI / 2)
    end
    lp
end

"Return the ctsem log posterior and its raw-coordinate gradient."
function ctsem_logposterior_and_gradient(rawest_r, sd; chunksize::Integer=7,
                                         threaded::Bool=false)
    x = fvec(rawest_r)
    nparams = Int(field(sd, :nparams))
    ntie = Int(field(sd, :ntipredeffects))
    priormod = field(sd, :priormod) / Int(field(sd, :nsubsets))
    scale = field(sd, :tipredeffectscale)
    lp = ctsem_logprior(x, sd)
    ll = ctsem_loglik_by_subject(x, sd; threaded=threaded)
    # ctsem_gradient_by_subject returns the negative-log-likelihood gradient.
    g = -ctsem_gradient_by_subject(
        x, sd; chunksize=chunksize, progress_every=0, threaded=threaded)
    g[1:nparams] .-= priormod .* @view(x[1:nparams])
    ntie > 0 && (g[nparams+1:end] .-=
        (priormod / scale^2) .* @view(x[nparams+1:end]))
    (ll + lp, g)
end

"""
    ctsem_nuts(rawest, sd; samples=1000, warmup=1000, chains=4, ...)

Run independent AdvancedHMC NUTS chains with Stan-style diagonal or dense
metric and step-size adaptation. Returns a `samples × parameters × chains`
array and per-draw sampler diagnostics. R computes rank-normalized R-hat and
ESS.
"""
function ctsem_nuts(rawest_r, sd; samples::Integer=1000,
                    warmup::Integer=1000, chains::Integer=4,
                    target_accept::Real=0.8, max_depth::Integer=10,
                    chunksize::Integer=7, threaded::Bool=false,
                    seed::Integer=1, init_jitter::Real=0.01,
                    progress::Bool=true, metric::Symbol=:diagonal,
                    study_re::Bool=false, verbose::Bool=true,
                    progress_callback=nothing)
    samples > 0 || error("samples must be positive")
    warmup > 0 || error("warmup must be positive")
    chains > 0 || error("chains must be positive")
    0 < target_accept < 1 || error("target_accept must lie between zero and one")
    max_depth > 0 || error("max_depth must be positive")
    init_jitter >= 0 || error("init_jitter must be non-negative")
    metric in (:diagonal, :dense) ||
        error("metric must be :diagonal or :dense")
    x0 = Float64.(fvec(rawest_r))
    # Validate the prior branch before starting expensive adaptation.
    if study_re
        ctsem_re_logposterior_and_gradient(
            x0, sd; chunksize=chunksize, threaded=false)
    else
        ctsem_logprior(x0, sd)
    end
    D = length(x0)
    draws = Array{Float64}(undef, samples, D, chains)
    acceptance_rate = Matrix{Float64}(undef, samples, chains)
    tree_depth = Matrix{Int}(undef, samples, chains)
    n_steps = Matrix{Int}(undef, samples, chains)
    divergent = Matrix{Bool}(undef, samples, chains)
    warmup_divergent = Matrix{Bool}(undef, warmup, chains)
    step_size = Matrix{Float64}(undef, samples, chains)

    for chain in 1:chains
        verbose && println("NUTS chain ", chain, "/", chains)
        rng = Random.MersenneTwister(seed + chain - 1)
        initial = x0 + init_jitter .* randn(rng, D)
        posterior_and_gradient = study_re ?
            (x -> ctsem_re_logposterior_and_gradient(
                x, sd; chunksize=chunksize, threaded=threaded)) :
            (x -> ctsem_logposterior_and_gradient(
                x, sd; chunksize=chunksize, threaded=threaded))
        logdensity = x -> posterior_and_gradient(x)[1]
        logdensity_gradient = posterior_and_gradient
        metric_object = metric === :dense ?
            AdvancedHMC.DenseEuclideanMetric(D) :
            AdvancedHMC.DiagEuclideanMetric(D)
        hamiltonian = AdvancedHMC.Hamiltonian(
            metric_object, logdensity, logdensity_gradient)
        epsilon = AdvancedHMC.find_good_stepsize(rng, hamiltonian, initial)
        integrator = AdvancedHMC.Leapfrog(epsilon)
        kernel = AdvancedHMC.HMCKernel(
            AdvancedHMC.Trajectory{AdvancedHMC.MultinomialTS}(
                integrator,
                AdvancedHMC.GeneralisedNoUTurn(max_depth, 1000.0)))
        adaptor = AdvancedHMC.StanHMCAdaptor(
            AdvancedHMC.MassMatrixAdaptor(metric_object),
            AdvancedHMC.StepSizeAdaptor(target_accept, integrator))
        progress_enabled = progress || progress_callback !== nothing
        progress_update! = if progress_callback === nothing
            AdvancedHMC.simple_pm_next!
        else
            function (pm, stat)
                progress_callback(stat)
                progress && AdvancedHMC.simple_pm_next!(pm, stat)
                nothing
            end
        end
        chain_draws, stats = AdvancedHMC.sample(
            rng, hamiltonian, kernel, initial, warmup + samples,
            adaptor, warmup; drop_warmup=false, verbose=verbose,
            progress=progress_enabled, pm_next! = progress_update!)
        for iteration in 1:warmup
            warmup_divergent[iteration, chain] = stats[iteration].numerical_error
        end
        for iteration in 1:samples
            retained = warmup + iteration
            draws[iteration, :, chain] .= chain_draws[retained]
            acceptance_rate[iteration, chain] = stats[retained].acceptance_rate
            tree_depth[iteration, chain] = stats[retained].tree_depth
            n_steps[iteration, chain] = stats[retained].n_steps
            divergent[iteration, chain] = stats[retained].numerical_error
            step_size[iteration, chain] = stats[retained].step_size
        end
    end
    (
        draws=draws,
        acceptance_rate=acceptance_rate,
        tree_depth=tree_depth,
        n_steps=n_steps,
        divergent=divergent,
        warmup_divergent=warmup_divergent,
        step_size=step_size,
        target_accept=Float64(target_accept),
        max_depth=Int(max_depth),
        warmup=Int(warmup),
        samples=Int(samples),
        chains=Int(chains),
        seed=Int(seed),
        metric=metric,
        study_re=study_re
    )
end

"""
    ctsem_nuts_parallel(rawest, sd; chains=4, parallel_chains=2,
                        threads_per_chain=1, source_file=@__FILE__, ...)

Run independent NUTS chains in separate local Julia worker processes. Each
worker has its own Julia thread pool, random-number generator, metric adaptor,
and AdvancedHMC state. Results are returned in the same chain order and format
as `ctsem_nuts`, so fixed seeds reproduce the sequential-chain results.
"""
function ctsem_nuts_parallel(rawest_r, sd; samples::Integer=1000,
                             warmup::Integer=1000, chains::Integer=4,
                             parallel_chains::Integer=2,
                             threads_per_chain::Integer=1,
                             blas_threads::Integer=1,
                             target_accept::Real=0.8,
                             max_depth::Integer=10,
                             chunksize::Integer=7,
                             seed::Integer=1,
                             init_jitter::Real=0.01,
                             metric::Symbol=:diagonal,
                             study_re::Bool=false,
                             source_file::AbstractString=@__FILE__,
                             project_dir::AbstractString=@__DIR__,
                             progress::Bool=true,
                             progress_every::Integer=25)
    chains > 0 || error("chains must be positive")
    parallel_chains > 0 || error("parallel_chains must be positive")
    threads_per_chain > 0 || error("threads_per_chain must be positive")
    blas_threads > 0 || error("blas_threads must be positive")
    progress_every > 0 || error("progress_every must be positive")
    workers_requested = min(Int(parallel_chains), Int(chains))
    if workers_requested == 1
        return ctsem_nuts(
            rawest_r, sd; samples=samples, warmup=warmup, chains=chains,
            target_accept=target_accept, max_depth=max_depth,
            chunksize=chunksize, threaded=threads_per_chain > 1,
            seed=seed, init_jitter=init_jitter, progress=progress,
            metric=metric, study_re=study_re)
    end

    source_path = abspath(source_file)
    project_path = abspath(project_dir)
    isfile(source_path) || error("Cannot find Julia backend at $source_path")
    x0 = Float64.(fvec(rawest_r))
    # JuliaCall represents R data.frames as DataFrames.DataFrame.  Convert
    # these to ordinary matrices before worker serialization so workers need
    # only the backend project rather than JuliaCall's auxiliary packages.
    portable_sd = Dict{Symbol,Any}()
    for (key, value) in pairs(sd)
        portable_sd[Symbol(key)] =
            string(parentmodule(typeof(value))) == "DataFrames" ?
            Matrix(value) : value
    end
    flags = `--project=$(project_path) --threads=$(Int(threads_per_chain))`
    worker_ids = Distributed.addprocs(workers_requested; exeflags=flags)
    started = time()
    println(
        "Running ", chains, " NUTS chains on ", workers_requested,
        " Julia worker processes (", threads_per_chain,
        " thread(s) per chain)."
    )
    try
        for worker in worker_ids
            Distributed.remotecall_fetch(
                Core.eval, worker, Main,
                :(Base.include(Main, $source_path))
            )
            Distributed.remotecall_fetch(
                Core.eval, worker, Main, :(ctma_parallel_x0 = $x0)
            )
            Distributed.remotecall_fetch(
                Core.eval, worker, Main, :(ctma_parallel_sd = $portable_sd)
            )
            Distributed.remotecall_fetch(
                Core.eval, worker, Main,
                :(CtsemJulia.LinearAlgebra.BLAS.set_num_threads($blas_threads))
            )
        end
        available_workers = Channel{Int}(workers_requested)
        foreach(worker -> put!(available_workers, worker), worker_ids)
        total_iterations = Int(warmup + samples)
        progress_channel = progress ? Distributed.RemoteChannel(
            () -> Channel{Any}(max(64, Int(chains) * 8))
        ) : nothing
        monitor = nothing
        if progress
            monitor = @async begin
                latest = zeros(Int, Int(chains))
                completed = falses(Int(chains))
                completed_chains = 0
                while completed_chains < chains
                    message = take!(progress_channel)
                    chain = Int(message.chain)
                    if message.kind === :iteration
                        latest[chain] = max(
                            latest[chain], Int(message.iteration))
                    elseif message.kind === :done
                        latest[chain] = total_iterations
                        if !completed[chain]
                            completed[chain] = true
                            completed_chains += 1
                        end
                    end
                    achieved = sum(latest)
                    target = Int(chains) * total_iterations
                    elapsed = max(time() - started, eps(Float64))
                    rate = achieved / elapsed
                    eta = rate > 0 ? max(
                        (target - achieved) / rate, 0.0) : Inf
                    eta_text = isfinite(eta) ?
                        string(round(eta; digits=0), " s") : "calculating"
                    print(
                        "\rNUTS progress: ", achieved, "/", target,
                        " iterations (",
                        round(100 * achieved / target; digits=1), "%); ",
                        completed_chains, "/", chains,
                        " chains complete; elapsed ",
                        round(elapsed; digits=0), " s; ETA ", eta_text, "   "
                    )
                    flush(stdout)
                end
                println()
            end
        end
        tasks = map(collect(1:Int(chains))) do chain
            @async begin
                worker = take!(available_workers)
                try
                    callback = progress ? (stat -> begin
                        iteration = Int(stat.iterations)
                        if iteration == 1 || iteration == total_iterations ||
                           iteration % Int(progress_every) == 0
                            put!(progress_channel, (
                                kind=:iteration, chain=chain,
                                iteration=iteration
                            ))
                        end
                        nothing
                    end) : nothing
                    result = Distributed.remotecall_fetch(
                        Core.eval, worker, Main,
                        :(CtsemJulia.ctsem_nuts(
                            ctma_parallel_x0, ctma_parallel_sd;
                            samples=$samples, warmup=$warmup, chains=1,
                            target_accept=$target_accept, max_depth=$max_depth,
                            chunksize=$chunksize,
                            threaded=$(threads_per_chain > 1),
                            seed=$(seed + chain - 1), init_jitter=$init_jitter,
                            progress=false, metric=$(QuoteNode(metric)),
                            study_re=$study_re, verbose=false,
                            progress_callback=$callback))
                    )
                    result
                finally
                    progress && put!(progress_channel, (
                        kind=:done, chain=chain, iteration=total_iterations
                    ))
                    put!(available_workers, worker)
                end
            end
        end
        results = fetch.(tasks)
        progress && wait(monitor)
        progress && println(
            "Finished ", chains, " parallel NUTS chains in ",
            round(time() - started; digits=1), " s."
        )
        return (
            draws=cat((result.draws for result in results)...; dims=3),
            acceptance_rate=hcat((result.acceptance_rate for result in results)...),
            tree_depth=hcat((result.tree_depth for result in results)...),
            n_steps=hcat((result.n_steps for result in results)...),
            divergent=hcat((result.divergent for result in results)...),
            warmup_divergent=hcat(
                (result.warmup_divergent for result in results)...),
            step_size=hcat((result.step_size for result in results)...),
            target_accept=Float64(target_accept),
            max_depth=Int(max_depth),
            warmup=Int(warmup),
            samples=Int(samples),
            chains=Int(chains),
            seed=Int(seed),
            metric=metric,
            study_re=study_re,
            parallel_chains=workers_requested,
            threads_per_chain=Int(threads_per_chain),
            blas_threads=Int(blas_threads)
        )
    finally
        Distributed.rmprocs(worker_ids)
    end
end

"""
    ctsem_optimize_multistart(rawest, sd; starts=4, parallel_starts=2,
                              threads_per_start=1, ...)

Run independent automatically scaled L-BFGS optimizations in separate local
Julia worker processes. Start 1 uses the supplied vector exactly; later starts
receive reproducible relative raw-scale jitter. The best converged solution is
selected by log likelihood, with the best finite solution retained when none
meets the requested raw-gradient tolerance.
"""
function ctsem_optimize_multistart(rawest_r, sd; starts::Integer=4,
                                   parallel_starts::Integer=2,
                                   threads_per_start::Integer=1,
                                   blas_threads::Integer=1,
                                   start_jitter::Real=0.1,
                                   seed::Integer=1,
                                   staged::Bool=true,
                                   maxiters::Integer=100,
                                   chunksize::Integer=7,
                                   initial_raw_step::Real=0.01,
                                   gradient_tolerance::Real=1e-3,
                                   function_relative_tolerance::Real=0.0,
                                   parameter_relative_tolerance::Real=0.0,
                                   source_file::AbstractString=@__FILE__,
                                   project_dir::AbstractString=@__DIR__,
                                   progress::Bool=true)
    starts > 0 || error("starts must be positive")
    parallel_starts > 0 || error("parallel_starts must be positive")
    threads_per_start > 0 || error("threads_per_start must be positive")
    blas_threads > 0 || error("blas_threads must be positive")
    start_jitter >= 0 || error("start_jitter must be non-negative")
    workers_requested = min(Int(parallel_starts), Int(starts))
    x0 = Float64.(fvec(rawest_r))
    jitter_scale = max.(abs.(x0), 0.1)
    start_vectors = Vector{Vector{Float64}}(undef, starts)
    start_vectors[1] = copy(x0)
    for start in 2:starts
        rng = Random.MersenneTwister(seed + start - 1)
        start_vectors[start] = x0 .+
            Float64(start_jitter) .* jitter_scale .* randn(rng, length(x0))
    end

    source_path = abspath(source_file)
    project_path = abspath(project_dir)
    isfile(source_path) || error("Cannot find Julia backend at $source_path")
    portable_sd = Dict{Symbol,Any}()
    for (key, value) in pairs(sd)
        portable_sd[Symbol(key)] =
            string(parentmodule(typeof(value))) == "DataFrames" ?
            Matrix(value) : value
    end
    flags = `--project=$(project_path) --threads=$(Int(threads_per_start))`
    worker_ids = Distributed.addprocs(workers_requested; exeflags=flags)
    started = time()
    println(
        "Running ", starts, " ML starts on ", workers_requested,
        " Julia worker processes (", threads_per_start,
        " thread(s) per start)."
    )
    try
        for worker in worker_ids
            Distributed.remotecall_fetch(
                Core.eval, worker, Main, :(Base.include(Main, $source_path)))
            Distributed.remotecall_fetch(
                Core.eval, worker, Main, :(ctma_parallel_sd = $portable_sd))
            Distributed.remotecall_fetch(
                Core.eval, worker, Main,
                :(CtsemJulia.LinearAlgebra.BLAS.set_num_threads($blas_threads)))
        end
        available_workers = Channel{Int}(workers_requested)
        foreach(worker -> put!(available_workers, worker), worker_ids)
        tasks = map(collect(1:Int(starts))) do start
            @async begin
                worker = take!(available_workers)
                try
                    start_time = time()
                    start_vector = start_vectors[start]
                    result = Distributed.remotecall_fetch(
                        Core.eval, worker, Main,
                        quote
                            local xstart = $start_vector
                            try
                                local opt = if $staged
                                    CtsemJulia.ctsem_optimize_staged(
                                        xstart, ctma_parallel_sd;
                                        maxiters=$maxiters,
                                        chunksize=$chunksize,
                                        show_trace=false,
                                        threaded=$(threads_per_start > 1),
                                        initial_raw_step=$initial_raw_step,
                                        gradient_tolerance=$gradient_tolerance,
                                        function_relative_tolerance=
                                            $function_relative_tolerance,
                                        parameter_relative_tolerance=
                                            $parameter_relative_tolerance)
                                else
                                    CtsemJulia.ctsem_optimize_auto(
                                        xstart, ctma_parallel_sd;
                                        maxiters=$maxiters,
                                        chunksize=$chunksize,
                                        show_trace=false,
                                        threaded=$(threads_per_start > 1),
                                        initial_raw_step=$initial_raw_step,
                                        gradient_tolerance=$gradient_tolerance,
                                        function_relative_tolerance=
                                            $function_relative_tolerance,
                                        parameter_relative_tolerance=
                                            $parameter_relative_tolerance)
                                end
                                local x = opt.estimates
                                local g = CtsemJulia.ctsem_gradient_by_subject(
                                    x, ctma_parallel_sd;
                                    chunksize=$chunksize, progress_every=0,
                                    threaded=$(threads_per_start > 1))
                                local maxg = maximum(abs.(g))
                                (
                                    success=true, error="", estimates=x,
                                    loglik=-CtsemJulia.Optim.minimum(opt.result),
                                    iterations=opt.total_iterations,
                                    stage_iterations=opt.stage_iterations,
                                    stage_objectives=opt.stage_objectives,
                                    stage_maximum_absolute_gradients=
                                        opt.stage_maximum_absolute_gradients,
                                    optimizer_converged=
                                        CtsemJulia.Optim.converged(opt.result),
                                    termination_code=string(
                                        opt.result.termination_code),
                                    stopped_by=string(opt.result.stopped_by),
                                    converged=maxg <= $gradient_tolerance,
                                    initial_maximum_absolute_gradient=
                                        opt.initial_maximum_absolute_gradient,
                                    maximum_absolute_gradient=maxg,
                                    gradient_L2_norm=sqrt(sum(g .^ 2)),
                                    maximum_parameter_change=
                                        maximum(abs.(x .- xstart)),
                                    parameter_change_L2=
                                        sqrt(sum((x .- xstart) .^ 2)),
                                    scale_minimum=minimum(opt.scale),
                                    scale_maximum=maximum(opt.scale),
                                    start_vector=xstart
                                )
                            catch err
                                (
                                    success=false,
                                    error=sprint(showerror, err),
                                    estimates=fill(NaN, length(xstart)),
                                    loglik=-Inf, iterations=0,
                                    stage_iterations=Int[],
                                    stage_objectives=Float64[],
                                    stage_maximum_absolute_gradients=Float64[],
                                    optimizer_converged=false,
                                    termination_code="Error",
                                    stopped_by="Error", converged=false,
                                    initial_maximum_absolute_gradient=Inf,
                                    maximum_absolute_gradient=Inf,
                                    gradient_L2_norm=Inf,
                                    maximum_parameter_change=Inf,
                                    parameter_change_L2=Inf,
                                    scale_minimum=NaN, scale_maximum=NaN,
                                    start_vector=xstart
                                )
                            end
                        end
                    )
                    progress && println(
                        "ML start ", start, "/", starts, " finished in ",
                        round(time() - start_time; digits=1), " s; logLik ",
                        round(result.loglik; digits=6), "; max gradient ",
                        round(result.maximum_absolute_gradient; sigdigits=5),
                        result.converged ? " (converged)." : " (not converged)."
                    )
                    result
                finally
                    put!(available_workers, worker)
                end
            end
        end
        results = fetch.(tasks)
        successful = findall(result -> result.success && isfinite(result.loglik),
                             results)
        isempty(successful) && error(
            "All Julia ML starts failed: " *
            join(unique(result.error for result in results), " | "))
        converged = filter(start -> results[start].converged, successful)
        candidates = isempty(converged) ? successful : converged
        best_index = candidates[argmax(results[start].loglik for start in candidates)]
        progress && println(
            "Selected ML start ", best_index, "/", starts, " with logLik ",
            round(results[best_index].loglik; digits=8), "."
        )
        diagnostics = (
            start=collect(1:Int(starts)),
            success=[result.success for result in results],
            converged=[result.converged for result in results],
            loglik=[result.loglik for result in results],
            maximum_absolute_gradient=[
                result.maximum_absolute_gradient for result in results],
            iterations=[result.iterations for result in results],
            termination_code=[result.termination_code for result in results],
            error=[result.error for result in results]
        )
        return (
            best=results[best_index],
            best_start=best_index,
            diagnostics=diagnostics,
            starts=Int(starts),
            parallel_starts=workers_requested,
            threads_per_start=Int(threads_per_start),
            start_jitter=Float64(start_jitter),
            seed=Int(seed),
            elapsed=time() - started
        )
    finally
        Distributed.rmprocs(worker_ids)
    end
end

"""
    ctsem_hessian_by_subject(rawest, sd; chunksize=7,
                             progress_every=1000, threaded=false)

Compute the observed Hessian of the negative log-likelihood. Each subject is
differentiated on the much smaller local population-parameter scale. Its local
Hessian is then mapped exactly to the global population and TI-effect parameter
vector. Because that map is linear, no second-derivative mapping term is needed.
"""
function ctsem_hessian_by_subject(rawest_r, sd; chunksize::Integer=7,
                                  progress_every::Integer=1000,
                                  threaded::Bool=false)
    check_supported(sd)
    x = fvec(rawest_r)
    d = prepare_subject_data(sd)
    rawpop = x[1:d.nparams]
    tipar = x[d.nparams+1:end]
    TIP = zeros(eltype(x), d.nparams, d.ntipred)
    for j in axes(TIP, 2), p in axes(TIP, 1)
        ix = d.tie_setup[p, j]
        ix > 0 && (TIP[p, j] = tipar[ix])
    end

    global_indices = Vector{Vector{Int}}(undef, d.nparams)
    predictor_indices = Vector{Vector{Int}}(undef, d.nparams)
    for p in 1:d.nparams
        gi = Int[p]
        pj = Int[0]
        for j in 1:d.ntipred
            ix = d.tie_setup[p, j]
            if ix > 0
                push!(gi, d.nparams + ix)
                push!(pj, j)
            end
        end
        global_indices[p] = gi
        predictor_indices[p] = pj
    end

    function add_subject_hessian!(H, si)
        z = vec(@view d.tipreds[si, :])
        localraw = rawpop + TIP * z
        f = q -> -subject_loglik(q, d, d.firstrow[si], d.lastrow[si])
        local_chunksize = min(chunksize, length(localraw))
        cfg = ForwardDiff.HessianConfig(
            f, localraw, ForwardDiff.Chunk{local_chunksize}())
        hs = ForwardDiff.hessian(f, localraw, cfg)

        for p in 1:d.nparams, q in 1:d.nparams
            hpq = hs[p, q]
            iszero(hpq) && continue
            gip = global_indices[p]
            gjq = global_indices[q]
            pjp = predictor_indices[p]
            pjq = predictor_indices[q]
            for a in eachindex(gip)
                ca = pjp[a] == 0 ? one(eltype(x)) : z[pjp[a]]
                ia = gip[a]
                for b in eachindex(gjq)
                    cb = pjq[b] == 0 ? one(eltype(x)) : z[pjq[b]]
                    H[ia, gjq[b]] += hpq * ca * cb
                end
            end
        end
        nothing
    end

    if threaded && Threads.nthreads() > 1
        thread_H = [zeros(eltype(x), length(x), length(x))
                    for _ in 1:Threads.maxthreadid()]
        completed = Threads.Atomic{Int}(0)
        Threads.@threads :static for si in 1:d.nsubjects
            add_subject_hessian!(thread_H[Threads.threadid()], si)
            done = Threads.atomic_add!(completed, 1) + 1
            progress_every > 0 && done % progress_every == 0 &&
                println("subject-Hessian progress: ", done, "/", d.nsubjects)
        end
        H = zeros(eltype(x), length(x), length(x))
        for Ht in thread_H
            H .+= Ht
        end
    else
        H = zeros(eltype(x), length(x), length(x))
        for si in 1:d.nsubjects
            add_subject_hessian!(H, si)
            progress_every > 0 && si % progress_every == 0 &&
                println("subject-Hessian progress: ", si, "/", d.nsubjects)
        end
    end
    Matrix(Symmetric((H + H') / 2))
end

"Compute the observed-information covariance and curvature diagnostics."
function ctsem_vcov_by_subject(rawest, sd; chunksize::Integer=7,
                               progress_every::Integer=1000,
                               threaded::Bool=false,
                               absolute_eigenvalue_tolerance::Real=1e-8,
                               relative_eigenvalue_tolerance::Real=1e-12)
    H = ctsem_hessian_by_subject(
        rawest, sd; chunksize=chunksize,
        progress_every=progress_every, threaded=threaded)
    eigenvalues = eigvals(Symmetric(H))
    minimum_eigenvalue = minimum(eigenvalues)
    maximum_eigenvalue = maximum(eigenvalues)
    number_nonpositive = count(<=(0), eigenvalues)
    eigenvalue_tolerance = max(
        Float64(absolute_eigenvalue_tolerance),
        Float64(relative_eigenvalue_tolerance) * maximum(abs, eigenvalues))
    minimum_eigenvalue >= -eigenvalue_tolerance || error(
        "Observed Hessian is not positive definite: " *
        "$(number_nonpositive) non-positive eigenvalue(s); " *
        "minimum=$(minimum_eigenvalue), numerical tolerance=" *
        "$(eigenvalue_tolerance)")

    # Subject-wise accumulation can leave an otherwise zero-curvature direction
    # a few ulps below zero.  In that narrow numerical boundary region, use a
    # symmetric eigendecomposition and floor the affected eigenvalues before
    # inversion.  Materially negative curvature is rejected above.
    decomposition = eigen(Symmetric(H))
    regularized_eigenvalues = max.(decomposition.values, eigenvalue_tolerance)
    number_regularized = count(<(eigenvalue_tolerance), decomposition.values)
    C = decomposition.vectors * Diagonal(inv.(regularized_eigenvalues)) *
        decomposition.vectors'
    (
        hessian=H,
        covariance=Matrix(Symmetric(C)),
        minimum_eigenvalue=minimum_eigenvalue,
        maximum_eigenvalue=maximum_eigenvalue,
        condition_number=maximum_eigenvalue / minimum(regularized_eigenvalues),
        number_nonpositive=number_nonpositive,
        number_regularized=number_regularized,
        eigenvalue_tolerance=eigenvalue_tolerance,
        regularized_minimum_eigenvalue=minimum(regularized_eigenvalues)
    )
end

"L-BFGS using the validated subject-wise score decomposition."
function ctsem_optimize_by_subject(rawest, sd; maxiters::Integer=5,
                                   chunksize::Integer=7,
                                   show_trace::Bool=true,
                                   threaded::Bool=false)
    x = fvec(rawest)
    function f(z)
        try
            value = if threaded
                -ctsem_loglik_by_subject(z, sd; threaded=true)
            else
                ctsem_negloglik(z, sd)
            end
            isfinite(value) ? value : Inf
        catch err
            # A threaded likelihood wraps a subject-level factorization error
            # in TaskFailedException. Such trial points are outside the valid
            # covariance region and should be rejected by the line search.
            if err isa PosDefException || err isa SingularException ||
               err isa DomainError || err isa TaskFailedException
                return Inf
            end
            rethrow()
        end
    end
    function g!(g, z)
        g .= ctsem_gradient_by_subject(z, sd; chunksize=chunksize,
                                      progress_every=0, threaded=threaded)
        nothing
    end
    options = Optim.Options(iterations=maxiters, show_trace=show_trace,
                            store_trace=true, allow_f_increases=false)
    # Raw ctsem coordinates can have very different curvature. In particular,
    # the neutral zero initialization has a large score, so an unscaled unit
    # trial step immediately leaves the positive-definite covariance region.
    initial_step = Optim.LineSearches.InitialStatic(alpha=1e-6, scaled=false)
    method = Optim.LBFGS(
        alphaguess=initial_step,
        linesearch=Optim.LineSearches.BackTracking(iterations=40))
    Optim.optimize(f, g!, x, method, options)
end

"""
    ctsem_optimize_auto(rawest, sd; ...)

Automatically precondition raw ctsem coordinates from the initial subject
score. The scale is chosen so that the first quasi-Newton direction proposes a
bounded raw displacement per parameter before backtracking. This is especially
important for neutral zero initialization, where variance scores can be many
orders of magnitude larger than other scores.
"""
function ctsem_optimize_auto(rawest, sd; maxiters::Integer=100,
                             chunksize::Integer=7,
                             show_trace::Bool=true,
                             threaded::Bool=false,
                             initial_raw_step::Real=0.01,
                             gradient_tolerance::Real=1e-3,
                             function_relative_tolerance::Real=0.0,
                             parameter_relative_tolerance::Real=0.0)
    x0 = fvec(rawest)
    initial_raw_step > 0 || error("initial_raw_step must be positive")
    gradient_tolerance > 0 || error("gradient_tolerance must be positive")
    function_relative_tolerance >= 0 ||
        error("function_relative_tolerance must be non-negative")
    parameter_relative_tolerance >= 0 ||
        error("parameter_relative_tolerance must be non-negative")

    g0 = ctsem_gradient_by_subject(
        x0, sd; chunksize=chunksize, progress_every=0, threaded=threaded)
    # With x = x0 + scale*u and the first u direction -scale*g, the raw
    # displacement is -scale^2*g. Bound that displacement by
    # `initial_raw_step`, while leaving near-flat coordinates on their natural
    # raw scale.
    scale = sqrt.(initial_raw_step ./ max.(abs.(g0), initial_raw_step))
    scale = clamp.(scale, 1e-6, 1.0)
    u0 = zeros(eltype(x0), length(x0))
    to_raw = u -> x0 + scale .* u

    function f(u)
        try
            value = if threaded
                -ctsem_loglik_by_subject(to_raw(u), sd; threaded=true)
            else
                ctsem_negloglik(to_raw(u), sd)
            end
            isfinite(value) ? value : Inf
        catch err
            if err isa PosDefException || err isa SingularException ||
               err isa DomainError || err isa TaskFailedException
                return Inf
            end
            rethrow()
        end
    end
    function g!(g, u)
        gx = ctsem_gradient_by_subject(
            to_raw(u), sd; chunksize=chunksize,
            progress_every=0, threaded=threaded)
        g .= scale .* gx
        nothing
    end

    # A scaled gradient below min(scale)*raw tolerance guarantees that every
    # raw gradient component is below the requested tolerance. Function and
    # parameter tolerances may stop earlier, but the R adapter independently
    # checks the final raw gradient before declaring convergence.
    scaled_g_tolerance = max(eps(Float64), minimum(scale) * gradient_tolerance)
    options = Optim.Options(
        iterations=maxiters,
        show_trace=false,
        store_trace=false,
        callback=optimizer_progress_callback(show_trace),
        allow_f_increases=false,
        g_abstol=scaled_g_tolerance,
        f_reltol=function_relative_tolerance,
        x_reltol=parameter_relative_tolerance)
    method = Optim.LBFGS(
        linesearch=Optim.LineSearches.BackTracking(iterations=40))
    result = Optim.optimize(f, g!, u0, method, options)
    show_trace && println()
    estimates = to_raw(Optim.minimizer(result))
    final_gradient = ctsem_gradient_by_subject(
        estimates, sd; chunksize=chunksize,
        progress_every=0, threaded=threaded)
    (result=result, estimates=estimates, scale=scale,
     initial_gradient=g0, initial_maximum_absolute_gradient=maximum(abs.(g0)),
     total_iterations=Optim.iterations(result),
     stage_iterations=[Optim.iterations(result)],
     stage_objectives=[Optim.minimum(result)],
     stage_maximum_absolute_gradients=[maximum(abs.(final_gradient))])
end

"Optimize selected raw coordinates while holding all others fixed."
function ctsem_optimize_active(rawest, active, sd; maxiters::Integer=100,
                               chunksize::Integer=7,
                               show_trace::Bool=true,
                               threaded::Bool=false,
                               initial_raw_step::Real=0.01,
                               gradient_tolerance::Real=1e-3,
                               function_relative_tolerance::Real=0.0,
                               parameter_relative_tolerance::Real=0.0)
    x0 = Float64.(fvec(rawest))
    active = Int.(collect(active))
    isempty(active) && error("At least one active coordinate is required")
    g0 = ctsem_gradient_by_subject(
        x0, sd; chunksize=chunksize, progress_every=0, threaded=threaded)
    scale = sqrt.(initial_raw_step ./
        max.(abs.(g0[active]), initial_raw_step))
    scale = clamp.(scale, 1e-6, 1.0)
    u0 = zeros(length(active))
    function to_raw(u)
        x = copy(x0)
        x[active] .= x0[active] .+ scale .* u
        x
    end
    function f(u)
        try
            value = -ctsem_loglik_by_subject(
                to_raw(u), sd; threaded=threaded)
            isfinite(value) ? value : Inf
        catch err
            if err isa PosDefException || err isa SingularException ||
               err isa DomainError || err isa TaskFailedException
                return Inf
            end
            rethrow()
        end
    end
    function g!(g, u)
        gx = ctsem_gradient_by_subject(
            to_raw(u), sd; chunksize=chunksize,
            progress_every=0, threaded=threaded)
        g .= scale .* gx[active]
        nothing
    end
    options = Optim.Options(
        iterations=maxiters,
        show_trace=false,
        store_trace=false,
        callback=optimizer_progress_callback(show_trace),
        allow_f_increases=false,
        g_abstol=max(eps(Float64), minimum(scale) * gradient_tolerance),
        f_reltol=function_relative_tolerance,
        x_reltol=parameter_relative_tolerance)
    method = Optim.LBFGS(
        linesearch=Optim.LineSearches.BackTracking(iterations=40))
    result = Optim.optimize(f, g!, u0, method, options)
    show_trace && println()
    estimates = to_raw(Optim.minimizer(result))
    final_gradient = ctsem_gradient_by_subject(
        estimates, sd; chunksize=chunksize,
        progress_every=0, threaded=threaded)
    (result=result, estimates=estimates, scale=scale, active=active,
     initial_gradient=g0, final_gradient=final_gradient)
end

"""
    ctsem_optimize_staged(rawest, sd; ...)

For WEC/TI-predictor models, first optimize population parameters while the
predictor effects remain fixed, then introduce all effects, and finally restart
the joint optimization with scaling recomputed at the new point. `maxiters`
applies to each stage; a converged joint second stage skips the third stage.
"""
function ctsem_optimize_staged(rawest, sd; maxiters::Integer=100,
                               chunksize::Integer=7,
                               show_trace::Bool=true,
                               threaded::Bool=false,
                               initial_raw_step::Real=0.01,
                               gradient_tolerance::Real=1e-3,
                               function_relative_tolerance::Real=0.0,
                               parameter_relative_tolerance::Real=0.0)
    x0 = Float64.(fvec(rawest))
    nparams = Int(field(sd, :nparams))
    ntie = Int(field(sd, :ntipredeffects))
    ntie > 0 || return ctsem_optimize_auto(
        x0, sd; maxiters=maxiters, chunksize=chunksize,
        show_trace=show_trace, threaded=threaded,
        initial_raw_step=initial_raw_step,
        gradient_tolerance=gradient_tolerance,
        function_relative_tolerance=function_relative_tolerance,
        parameter_relative_tolerance=parameter_relative_tolerance)

    show_trace && println(
        "Optimization stage 1/3: population parameters; TI effects fixed")
    stages = Any[]
    stage1 = ctsem_optimize_active(
        x0, 1:nparams, sd; maxiters=maxiters, chunksize=chunksize,
        show_trace=show_trace, threaded=threaded,
        initial_raw_step=initial_raw_step,
        gradient_tolerance=gradient_tolerance,
        function_relative_tolerance=function_relative_tolerance,
        parameter_relative_tolerance=parameter_relative_tolerance)
    push!(stages, stage1)

    show_trace && println(
        "Optimization stage 2/3: population and TI-effect parameters")
    stage2 = ctsem_optimize_active(
        stage1.estimates, eachindex(x0), sd; maxiters=maxiters,
        chunksize=chunksize, show_trace=show_trace, threaded=threaded,
        initial_raw_step=initial_raw_step,
        gradient_tolerance=gradient_tolerance,
        function_relative_tolerance=function_relative_tolerance,
        parameter_relative_tolerance=parameter_relative_tolerance)
    push!(stages, stage2)
    final = stage2
    if maximum(abs.(stage2.final_gradient)) > gradient_tolerance
        show_trace && println(
            "Optimization stage 3/3: jointly rescaled refinement")
        final = ctsem_optimize_active(
            stage2.estimates, eachindex(x0), sd; maxiters=maxiters,
            chunksize=chunksize, show_trace=show_trace, threaded=threaded,
            initial_raw_step=initial_raw_step,
            gradient_tolerance=gradient_tolerance,
            function_relative_tolerance=function_relative_tolerance,
            parameter_relative_tolerance=parameter_relative_tolerance)
        push!(stages, final)
    end
    (
        result=final.result,
        estimates=final.estimates,
        scale=final.scale,
        initial_gradient=stage1.initial_gradient,
        initial_maximum_absolute_gradient=
            maximum(abs.(stage1.initial_gradient)),
        total_iterations=sum(Optim.iterations(s.result) for s in stages),
        stage_iterations=[Optim.iterations(s.result) for s in stages],
        stage_objectives=[Optim.minimum(s.result) for s in stages],
        stage_maximum_absolute_gradients=[
            maximum(abs.(s.final_gradient)) for s in stages]
    )
end

"Return raw-scale parameters from a scaled-coordinate Optim result."
scaled_minimizer(result, start, scale; scale_multiplier::Real=0.1) =
    fvec(start) + scale_multiplier .* fvec(scale) .* Optim.minimizer(result)

"""
    ctsem_optimize_scaled(start, sd, scale; ...)

Optimise `theta = start + scale .* u` from `u = 0`, using the validated
subject-wise raw-scale score and the chain rule `g_u = scale .* g_theta`.
`scale` is typically `sqrt.(diag(covariance))` from a prior fit.
"""
function ctsem_optimize_scaled(start, sd, scale; maxiters::Integer=5,
                               chunksize::Integer=7,
                               scale_multiplier::Real=0.1,
                               show_trace::Bool=true)
    x0 = fvec(start)
    s = scale_multiplier .* fvec(scale)
    length(s) == length(x0) || error("scale and start lengths differ")
    all(isfinite, s) || error("all scales must be finite")
    all(>(0), s) || error("all scales must be strictly positive")
    u0 = zeros(eltype(x0), length(x0))
    to_raw = u -> x0 + s .* u
    function f(u)
        try
            value = ctsem_negloglik(to_raw(u), sd)
            isfinite(value) ? value : Inf
        catch err
            if err isa PosDefException || err isa SingularException ||
               err isa DomainError
                return Inf
            end
            rethrow()
        end
    end
    function g!(g, u)
        gx = ctsem_gradient_by_subject(to_raw(u), sd;
                 chunksize=chunksize, progress_every=0)
        g .= s .* gx
        nothing
    end
    options = Optim.Options(iterations=maxiters, show_trace=show_trace,
                            store_trace=true, allow_f_increases=false)
    method = Optim.LBFGS(linesearch=Optim.LineSearches.BackTracking())
    Optim.optimize(f, g!, u0, method, options)
end

"Return raw parameters from a full-Cholesky-whitened Optim result."
whitened_minimizer(result, start, factor; scale_multiplier::Real=0.1) =
    fvec(start) + scale_multiplier .* factor * Optim.minimizer(result)

"""L-BFGS in full covariance-whitened coordinates `theta=start+m*L*u`."""
function ctsem_optimize_whitened(start, sd, factor; maxiters::Integer=5,
                                  chunksize::Integer=7,
                                  scale_multiplier::Real=0.1,
                                  show_trace::Bool=true)
    x0 = fvec(start)
    L = factor
    size(L) == (length(x0), length(x0)) ||
        error("factor dimensions must match start length")
    all(isfinite, L) || error("factor entries must be finite")
    u0 = zeros(eltype(x0), length(x0))
    to_raw = u -> x0 + scale_multiplier .* (L * u)
    function f(u)
        try
            value = ctsem_negloglik(to_raw(u), sd)
            isfinite(value) ? value : Inf
        catch err
            if err isa PosDefException || err isa SingularException ||
               err isa DomainError
                return Inf
            end
            rethrow()
        end
    end
    function g!(g, u)
        gx = ctsem_gradient_by_subject(to_raw(u), sd;
                 chunksize=chunksize, progress_every=0)
        g .= scale_multiplier .* (L' * gx)
        nothing
    end
    options = Optim.Options(iterations=maxiters, show_trace=show_trace,
                            store_trace=true, allow_f_increases=false)
    method = Optim.LBFGS(linesearch=Optim.LineSearches.BackTracking())
    Optim.optimize(f, g!, u0, method, options)
end

"Compute the ForwardDiff gradient of the negative log-likelihood."
function ctsem_gradient(rawest, sd; chunksize::Integer=8)
    x = fvec(rawest)
    f = z -> ctsem_negloglik(z, sd)
    cfg = ForwardDiff.GradientConfig(f, x, ForwardDiff.Chunk{chunksize}())
    ForwardDiff.gradient(f, x, cfg)
end

"Run L-BFGS with ForwardDiff gradients; intended only after gradient validation."
function ctsem_optimize(rawest, sd; maxiters::Integer=5,
                        chunksize::Integer=8, show_trace::Bool=true)
    x = fvec(rawest)
    f = z -> ctsem_negloglik(z, sd)
    cfg = ForwardDiff.GradientConfig(f, x, ForwardDiff.Chunk{chunksize}())
    function g!(g, z)
        ForwardDiff.gradient!(g, f, z, cfg)
        nothing
    end
    options = Optim.Options(iterations=maxiters, show_trace=show_trace,
                            store_trace=true, allow_f_increases=false)
    Optim.optimize(f, g!, x, Optim.LBFGS(), options)
end

end # module
