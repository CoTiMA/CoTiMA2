\name{ctmaCompute}
\alias{ctmaCompute}
\alias{ctmaHost}
\alias{ctmaCheckCompute}
\alias{ctmaJobStatus}
\alias{ctmaCollect}
\alias{ctmaResume}
\title{Experimental SSH execution plans and durable background jobs}
\description{Configure local and remote Julia workers, establish SSH tunnels, and submit or collect background fitting jobs. Configuration contains no passwords.}
\usage{
ctmaHost(host, user = NULL, workers = 1L, port = 22L,
  authentication = c("auto", "ask"), julia = "julia", project = NULL,
  directory = "~/CoTiMA2-Benchmark", rscript = "Rscript", depot = NULL,
  threads = NULL)
ctmaCompute(workers = c(local = 1L), threads = 1L, local = NULL,
  remote = NULL, tasks = c("auto", "likelihood", "starts", "studies"),
  tunnel = "ssh", reconnect = TRUE, checkpoint = TRUE,
  execution = c("foreground", "background"), coordinator = NULL,
  retries = 1L, directory = NULL)
ctmaCheckCompute(compute)
ctmaJobStatus(job)
ctmaCollect(job)
ctmaResume(job)
}
\arguments{
\item{host}{SSH alias, host name or IP address.}
\item{user}{Optional SSH user; otherwise use the SSH configuration.}
\item{workers}{Host worker count, or a uniquely named vector of counts for ctmaCompute. Zero local workers permits an external-only plan.}
\item{port}{SSH port.}
\item{authentication}{Use existing SSH keys/configuration (auto), or ask OpenSSH for interactive password authentication (ask). macOS uses a hidden-input OS dialogue; Linux requires a working terminal or SSH askpass configuration.}
\item{julia}{Remote Julia executable. Prefer an absolute path or a path starting with ~/ if remote noninteractive PATH lacks Julia.}
\item{project}{Remote Julia project directory containing the installed backend dependencies. NULL attempts discovery through the installed CoTiMA2 package.}
\item{directory}{Remote working directory for ctmaHost, or local cache/job directory for ctmaCompute. The remote directory must not contain spaces or shell metacharacters.}
\item{rscript}{Remote Rscript executable for package discovery and jobs.}
\item{depot}{Optional remote JULIA_DEPOT_PATH.}
\item{threads}{Threads per worker. A host-specific value overrides the plan default. BLAS stays at one thread. This is a resource budget, not CPU affinity.}
\item{local}{Optional override of the local worker count.}
\item{remote}{One ctmaHost or a list of host profiles.}
\item{tasks}{The current release implements starts and likelihood for Julia ML. auto selects starts, including within each sequential ctmaInit study. Use \code{"studies"} to distribute complete primary-study fits from \code{ctmaInit(backend = "ctsemjulia")}.}
\item{tunnel}{Currently only ssh.}
\item{reconnect}{Retry detected worker/connection losses. An expired password login may require a new interactive call; credentials are not retained.}
\item{checkpoint}{Persist validated optimization state. Restarting uses the best saved parameter vector and resets L-BFGS history, rather than promising bitwise continuation. Completed results can be reused only for matching data/code/settings.}
\item{execution}{Run immediately or return a persistent ctmaJob handle.}
\item{coordinator}{Optional remote ctmaHost owning the background R process. Remote background jobs require local=0 on the submitting machine. Workers naming the coordinator are converted to local workers there. Other hosts must be reachable/authenticatable from the coordinator.}
\item{retries}{Maximum additional worker-start/task attempts on detected network failure.}
\item{compute}{An execution plan.}
\item{job}{Saved ctmaJob handle.}
}
\details{
This experimental macOS/Linux implementation supports maximum likelihood without integrated subject/study random effects or priors. A shared likelihood is split over independent subjects; workers retain data between parameter evaluations. Multiple starts select a converged solution and fail explicitly when none converges. Unsupported settings are not silently approximated. With \code{tasks = "studies"}, each worker runs one primary-study fit at a time, with the thread budget specified for its host. Studies already loaded from disk are reused. Results are collected in the original study order and can be passed to \code{ctmaFit()}. This mode uses foreground R workers and reverse SSH tunnels. Install matching CoTiMA2 and ctsem builds on all hosts and use \code{rscript} to select the remote R environment. Background jobs, automatic reconnect and checkpoint recovery do not apply to this mode.

With \code{ctmaFit(backend = "ctsemjulia")}, \code{tasks = "likelihood"}
runs an experimental foreground ML optimization without priors. It supports
independent persons, including Laplace-integrated individual random effects.
Study-level random effects and missing time-independent predictors are outside
this route. Workers use matching Git-pinned ctsem installations and
noninteractive SSH. The coordinating R process optimizes the summed objective
and then passes the resulting raw parameter vector to the native ctsem optimizer
for full-data refinement and uncertainty. The final stage runs on the
coordinator using \code{coresToUse}. Worker optimization diagnostics are kept in
\code{distributed_likelihood}; the native fit supplies the final diagnostics.
This mode has no reconnect or checkpoint recovery.

OpenSSH builds encrypted tunnels. Host-key checking is retained. No password argument exists and R never receives the password. The ask dialogue is not suitable for unattended reconnection. Persistent remote background jobs use tmux; losing the submitting laptop does not stop the coordinator. Losing a required worker still requires retry/recovery, not merely a replacement tunnel.

Background jobs require options(cotima.source="/path/to/CoTiMA2") to create an exact Git source snapshot, plus Git, R and package dependencies on the coordinator. The snapshot is installed into a job-specific library without updating the normal package installation. Remote results, status and logs remain in the job directory. The handle is saved as job.rds in the local directory shown by print. Only owned worker processes and SSH masters are stopped; a background job is not stopped by status queries or disconnection.

For \code{backend = "ctsemjulia"}, \code{ctmaFit()} and
\code{ctmaInit()} accept a background plan with \code{local = 0},
\code{tasks = "auto"}, and one remote worker on the coordinator host.
The host thread setting supplies the native core budget. The coordinator
must already have a compatible ctsem and Julia environment. Each native
fit runs entirely on that host. Resuming a failed native job starts the
fit again.

Run ctmaCheckCompute before fitting to check SSH authentication and Julia availability. Worker startup additionally verifies Julia version and code hashes. This check does not promise that all models or credentials are valid. Resource controls in compute conflict with explicit backend parallel/thread controls; do not specify both. Source files are transferred via Git bundles, not by changing development checkouts.
}
\value{Host/plan configuration, a connectivity table, job status, a collected ordinary fitting result, or the resumed job handle.}
\examples{
plan <- ctmaCompute(workers = c(local = 0, studio1 = 4, studio2 = 4),
                    tasks = "likelihood")
\dontrun{
ctmaCheckCompute(plan)
fit <- ctmaFit(ctmaInitFit = init, backend = "julia", compute = plan)
status <- ctmaJobStatus(job)
result <- ctmaCollect(job)
}
}
