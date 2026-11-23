#' Create APA-style parameter tables
#'
#' @description Creates editable Word, legacy Word XML, PDF or LaTeX tables
#' from conventional CoTiMA and native ctsem matrix-parameter summaries.
#' It never fits a new model. Unsupported summary structures fail explicitly.
#' @param fit A fitted CoTiMA/ctsem object, or a named list of such objects
#'   for side-by-side comparisons. A \code{ctmaInit()} result produces rows
#'   for every primary study, using its retained aggregate summary.
#' @param summary Optional summary, or one summary per model. Its complete
#'   parameter table and explicit fit statistics must agree with the fit's
#'   stored or regenerated summary. This checks agreement, not the historical
#'   origin of an untagged legacy object. A fingerprint is attached on return.
#' @param parameters Character vector of exact parameter labels, matrix
#'   families, or case-insensitive selectors such as \code{"drift[2,1]"}.
#'   Bare matrix families such as \code{c(drift, diffusion)} are also accepted
#'   case-insensitively and always select every entry of those matrices.
#'   Character-vector variables remain supported; bare family names are reserved.
#'   Matrix rows are outcomes and columns are predictors. NULL asks yes/no
#'   for every available matrix entry in an interactive session.
#' @param statistics Character vector of source statistics or aliases, including
#'   \code{estimate}, \code{sd}, \code{median}, \code{"2.5"}, \code{"97.5"},
#'   \code{Rhat}, \code{ESS}, \code{bulk_ess}, and \code{tail_ess}. Only
#'   source-available statistics are used. Aggregate ctmaInit summaries do not
#'   retain medians or ESS/Rhat. SD is never silently relabelled SE;
#'   legacy ESS is not bulk ESS. NULL asks yes/no for each available statistic.
#' @param activeDirectory Existing output directory; defaults to the first
#'   fit's activeDirectory, otherwise the current working directory.
#' @param defaultOrder Order measurement, initial-state, mean, drift cross,
#'   drift auto, diffusion and remaining parameters. FALSE preserves selector order.
#' @param format One or more of \code{"docx"} (default), \code{"xml"}
#'   (Word 2003 WordprocessingML), \code{"tex"}, \code{"pdf"}; or
#'   \code{"none"} for an in-memory preview. Word requires optional packages
#'   officer/flextable; PDF requires XeLaTeX. Nothing is installed automatically.
#' @param file Safe output filename base, without a path or extension.
#' @param title Table title.
#' @param fontSize Positive font size in points for Word, XML, TeX and PDF output,
#'   default 10. Applies to table titles, headers, cells and notes. Plain data
#'   output has no font formatting.
#' @param digits Number of display decimals (0--12), default 4; numeric results retain source precision.
#' @param labels Optional named character vector replacing selected display labels;
#'   names identify parameter labels or matrix-coordinate selectors.
#' @param missing \code{"error"} for unavailable requested statistics, or
#'   \code{"note"} for blank cells with an explanatory note.
#' @param overwrite Replace existing requested outputs only when TRUE.
#' @param timeUnit Descriptive label for fitted time units; no automatic time conversion.
#' @param tableNote Additional character notes.
#' @param interactive Allow yes/no questions. Noninteractive and multiple-model
#'   calls require explicit parameters and statistics.
#' @param timeIntervals Optional distinct positive DT lags in fitted model time units.
#'   Select DRIFT entries to include their DT counterparts.
#' @param driftDraws Optional joint transformed drift array: draws by outcomes
#'   by predictors; one array per model for comparisons or per primary study
#'   for ctmaInit results. Cached Stan population DRIFT arrays and saved
#'   sampled Stan DRIFT draws are recognized automatically. Native Julia flat
#'   transformed draws are also recognized where applicable. User-supplied draws must
#'   correspond to the fitted matrix and are the user's responsibility.
#' @param dtPointOnly Explicitly permit DT plug-in point estimates without
#'   uncertainty. Otherwise joint draws or a regenerated native summary are required. CT interval endpoints
#'   are never exponentiated to manufacture DT intervals.
#' @param layout \code{"long"} (default) for portrait pages or
#'   \code{"wide"} for landscape pages. Independent of row arrangement.
#'   Columns that cannot fit legibly continue in labelled panels.
#' @param sortBy Case-insensitive row sorting: \code{"study"} (default, numeric
#'   study number), \code{"parameter"} (alphabetical parameter name), or a
#'   requested statistic/alias (ascending numeric order, missing values last).
#'   Every column moves with its row. An unrequested statistic is an error.
#'   Comparison tables require an unambiguous model-prefixed statistic column.
#'   Statistic sorting requires one parameter per row. Study sorting has no
#'   effect when there is no Study column. NULL preserves the original order.
#' @param sort Deprecated alias for \code{sortBy}. Do not supply both.
#' @details Use the parameter family `ricor` for individual random-intercept
#' correlations (CINT or MANIFESTMEANS). Add `rit0cor` for pairs involving
#' random T0MEANS. These selectors require a pooled fit or a single native fit
#' with an eligible saved correlation summary. Study-specific covariance-only
#' summaries cannot be used.
#' @param moderator Name or position of one moderator with coding stored by ctmaFit().
#' @param moderatorValues Optional moderator values in original units. Defaults to all categories or the mean and one standard deviation on either side. Study-level reference values give each study equal weight.
#' @param moderatorProfiles Optional matrix or data frame with one row per profile and columns named by the fitted TI predictors. Values use the fitted coding. Row names label the profiles. This also supports fits without stored moderator coding.
#' @param moderatorAt Named numeric vector of values for other TI predictors on the fitted scale. Omitted predictors are held at coded zero.
#' @param moderatorContrast Whether to report the second profile minus the first, computed on each joint draw. Requires exactly two profiles.
#' @param rowLayout \code{"parameter"} (default) puts one parameter and its
#'   statistics in each row. \code{"drift"} combines selected DRIFT entries
#'   into one row per study/model; DT drifts are grouped separately by lag.
#'   Other parameters retain separate rows. \code{"all"} puts all selected
#'   parameters in one row per study/model comparison. Compact columns carry
#'   parameter and statistic headers. Numeric data retain their long schema.
#' @param recomputeSummary Whether missing statistics or DT uncertainty may be
#'   obtained by calling summary() on saved native fits: \code{"ask"} (default),
#'   \code{"yes"} for batch use, or \code{"no"}. Uses
#'   \code{summary(saved_fit, timeinterval = lag)} and ctsem's model-specific
#'   transformations of retained joint draws. No model is refitted. A raw point
#'   estimate (typically stanfit$rawest) alone cannot provide uncertainty. Regular
#'   Stan optimization fits can reconstruct a missing transformed cache from
#'   rawposterior through their saved stanmodel and standata. No Julia is needed
#'   for regular Stan fits; no likelihood filtering or optimization is performed
#'   during this reconstruction.
#' @param modelNames Optional unique comparison-model labels, otherwise list names;
#'   for ctmaInit results, optional labels for all primary studies.
#' @details Table notes use explicit likelihood and parameter-count metadata;
#'   log-posterior values are never substituted for log-likelihood. MAP interval
#'   interpretation remains that of the stored approximation. Model comparisons
#'   align matrix coordinates and require matching substantive variable identities
#'   and time units. DT conversion assumes constant drift. It summarizes the
#'   chosen drift representation, not a mixture of individual trajectories.
#'   Notes distinguish sampling convergence from optimizer placement and report
#'   cached Laplace approximation status separately from optimizer convergence.
#'   Reading this metadata does not rerun an accuracy check. Stored
#'   native resultsSummary matrix tables can be used without a new summary call.
#'   Arbitrary hierarchical, empirical-moment, or raw random-effect summaries
#'   are not silently converted into matrix-parameter estimates. Obtain a
#'   supported transformed summary or select the relevant native fit first.
#' @return A \code{ctmaTables} list containing numeric \code{data}, formatted
#'   \code{display}, notes, output paths, validated summary, metadata, available
#'   matrix entries/statistics, and per-model tables for comparisons. The
#'   \code{layout} component records presentation options; \code{panels}
#'   contains the formatted column panels used for export.
#' @export
#' @section Moderator profiles:
#' A profile table shows conditional CT drift at each selected moderator value.
#' The calculation applies the model transformation after adding the TI effects
#' to the raw drift. Intervals and optional differences use the same joint raw
#' draws across profiles. Older fits and models with custom moderator transforms
#' or category contrasts can use explicit \code{moderatorProfiles}.
#'
#' @examples
#' \dontrun{
#' ctmaTables(initFit, parameters = c(drift, diffusion),
#'   statistics = c("estimate", "sd"), format = "none")
#' tables <- ctmaTables(fit, parameters = c("drift", "asymCINT"),
#'   statistics = c("estimate", "2.5", "median", "97.5"), format = "none")
#' ctmaTables(fit, moderator = "Burnout Measure", moderatorContrast = TRUE,
#'   format = "none")
#' ctmaTables(fit, parameters = "drift", statistics = c("estimate", "sd"),
#'   activeDirectory = tempdir(), format = c("docx", "tex"))
#' }

ctmaTables <- function(fit, summary=NULL, parameters=NULL, statistics=NULL,
    activeDirectory=NULL, defaultOrder=TRUE, format='docx', file='ctma_table',
    title='Model parameter estimates', digits=4, labels=NULL, missing=c('error','note'),
    overwrite=FALSE, timeUnit=NULL, tableNote=NULL, interactive=base::interactive(),
    timeIntervals=NULL, driftDraws=NULL, dtPointOnly=FALSE, modelNames=NULL,
    recomputeSummary=c("ask","yes","no"),
    layout=c("long","wide"), sortBy="study",
    rowLayout=c("parameter","drift","all"), fontSize=10, sort=NULL,
    moderator=NULL, moderatorValues=NULL, moderatorProfiles=NULL,
    moderatorAt=NULL, moderatorContrast=FALSE) {
  if (!base::missing(sort)) {
    if (!base::missing(sortBy)) .ctmaTablesStop('Supply only sortBy; sort is a deprecated alias.')
    .Deprecated(new='sortBy', package='CoTiMA2', old='sort')
    sortBy <- sort
  }
  parameters <- .ctmaTablesParameterInput(substitute(parameters),parent.frame())
  if(!is.numeric(fontSize)||length(fontSize)!=1L||!is.finite(fontSize)||fontSize<=0)
    .ctmaTablesStop("fontSize must be a single positive finite number of points.")
  if (!is.null(moderator) || !is.null(moderatorProfiles)) {
    if (!is.null(summary) || !is.null(timeIntervals))
      .ctmaTablesStop("Moderator profile tables require a fitted object and CT drift parameters.")
    if (!is.null(moderatorProfiles) && (!is.null(moderator) || !is.null(moderatorValues)))
      .ctmaTablesStop("Supply either moderator and optional moderatorValues or moderatorProfiles.")
    return(.ctmaTablesModerator(fit,moderator,moderatorValues,moderatorProfiles,
      moderatorAt,moderatorContrast,parameters,statistics,digits,format,
      activeDirectory,file,title,overwrite,tableNote,fontSize,defaultOrder))
  }
  if (!is.null(moderatorValues) || !is.null(moderatorAt) || !identical(moderatorContrast,FALSE))
    .ctmaTablesStop("Select moderator or moderatorProfiles for a moderator table.")
  missing <- match.arg(missing)
  recomputeSummary <- match.arg(recomputeSummary)
  layout <- match.arg(layout)
  rowLayout <- match.arg(rowLayout)
  known <- function(x) is.list(x)&&(inherits(x,c('CoTiMAFit','ctStanFit','ctJuliaFit','ctmaJuliaFit'))||is.list(x$summary))
  multi <- is.list(fit)&&!known(fit)&&length(fit)>0L&&all(vapply(fit,known,logical(1)))
  models <- if(multi)fit else list(fit)
  if(!is.logical(dtPointOnly)||length(dtPointOnly)!=1L||is.na(dtPointOnly)) .ctmaTablesStop('dtPointOnly must be TRUE or FALSE.')
  if(.ctmaTablesIsInit(fit)) {
    result <- .ctmaTablesInit(fit,summary,parameters,statistics,defaultOrder,digits,
      labels,missing,timeUnit,tableNote,interactive,timeIntervals,driftDraws,dtPointOnly,modelNames,recomputeSummary)
  } else {
  if(multi) {
    if(is.null(parameters)||is.null(statistics)) .ctmaTablesStop('Multiple-model tables require explicit parameters and statistics.')
    if(is.null(modelNames)) modelNames <- names(models)
    if(is.null(modelNames)) modelNames <- paste('Model',seq_along(models))
    if(!is.character(modelNames)||length(modelNames)!=length(models)||anyNA(modelNames)||any(!nzchar(modelNames))||anyDuplicated(modelNames)) .ctmaTablesStop('modelNames must uniquely label each fit.')
    if(!is.null(summary)&&(!is.list(summary)||length(summary)!=length(models))) .ctmaTablesStop('Supply one summary per model, in model order.')
    if(!is.null(driftDraws)&&(!is.list(driftDraws)||length(driftDraws)!=length(models))) .ctmaTablesStop('Supply one joint driftDraws array per model.')
  }
  tables <- lapply(seq_along(models),function(i) {
    supplied <- if(multi)summary[[i]] else summary
    current <- .ctmaTablesPrepare(models[[i]],supplied,parameters,statistics,timeIntervals,
      if(multi)driftDraws[[i]] else driftDraws,dtPointOnly,recomputeSummary,interactive,digits)
    r <- .ctmaTablesSingle(current,summary=if(is.null(current$ctma_tables_generated_summaries))supplied else NULL,
      parameters=parameters,statistics=statistics,activeDirectory=activeDirectory,defaultOrder=defaultOrder,
      format='none',file=file,title=title,digits=digits,labels=labels,missing=missing,
      overwrite=overwrite,timeUnit=timeUnit,tableNote=tableNote,interactive=interactive)
    if(!is.null(timeIntervals)) r <- .ctmaTablesAddDT(r,current,timeIntervals,if(multi)driftDraws[[i]] else driftDraws,dtPointOnly,digits,missing)
    r$generated_summaries <- current$ctma_tables_generated_summaries
    r
  })
  result <- tables[[1]]
  if(multi) {
    canonicalStats <- setdiff(names(tables[[1]]$data),c('matrix','row','col','parameter'))
    for(i in seq_along(tables)) {
      old <- setdiff(names(tables[[i]]$data),c('matrix','row','col','parameter'))
      if(length(old)!=length(canonicalStats)) .ctmaTablesStop('Model statistics do not have a consistent mapping.')
      names(tables[[i]]$data)[match(old,names(tables[[i]]$data))] <- canonicalStats
      names(tables[[i]]$display)[-1L] <- canonicalStats
    }
    ids <- lapply(tables,function(t)paste(t$data$matrix,t$data$row,t$data$col,sep=':'))
    allids <- unique(unlist(ids,use.names=FALSE))
    display <- data.frame(Parameter=rep('',length(allids)),stringsAsFactors=FALSE)
    long <- list();notes<-character()
    for(i in seq_along(tables)) {
      t<-tables[[i]];z<-match(ids[[i]],allids)
      blank <- !nzchar(display$Parameter[z]);display$Parameter[z[blank]]<-t$display$Parameter[blank]
      for(s in setdiff(names(t$display),'Parameter')) {
        v<-rep('',length(allids));v[z]<-t$display[[s]];display[[paste0(modelNames[i],': ',s)]]<-v
      }
      d<-t$data;d$Model<-modelNames[i];long[[i]]<-d
      notes<-c(notes,paste0(modelNames[i],': ',t$notes))
    }
    if(any(vapply(ids,length,integer(1))!=length(allids))) notes<-c(notes,'Blank model cells indicate a parameter absent from that model; no estimate was substituted.')
    result$data<-do.call(rbind,long);result$display<-display;result$notes<-notes
    result$metadata<-stats::setNames(lapply(tables,`[[`,'metadata'),modelNames)
    result$summary<-stats::setNames(lapply(tables,`[[`,'summary'),modelNames)
    result$tables<-stats::setNames(tables,modelNames)
    result$notes<-c(result$notes,'Models are aligned by matrix coordinates; comparison assumes the same variable identities and fitted time units in each model.')
  }
  }
  # Validate formats even if no writing is requested elsewhere in the workflow.
  if(!is.character(format)||!length(format)||anyNA(format)||any(!format %in% c('none','docx','xml','pdf','tex'))||('none' %in% format&&length(format)>1L)) .ctmaTablesStop('Unsupported output format.')
  # Presentation labels change only the display table, never numeric source columns.
  headings <- c(Mean='M',mean='M',sd='SD',Sd='SD',`50%`='Mdn',n_eff='ESS (legacy)',Rhat='Rhat')
  displayNames <- names(result$display)
  for(j in seq_along(displayNames)) for(s in names(headings)) {
    if(displayNames[j]==s) displayNames[j]<-headings[[s]]
    else if(endsWith(displayNames[j],paste0(': ',s))) displayNames[j]<-paste0(substr(displayNames[j],1,nchar(displayNames[j])-nchar(s)),headings[[s]])
  }
  names(result$display)<-displayNames
  if(any(grepl('Mdn',displayNames,fixed=TRUE))) result$notes<-c(result$notes,'Mdn = median.')
  if (!is.null(sortBy)) result <- .ctmaTablesSort(result, sortBy, multi, rowLayout)
  result <- .ctmaTablesLayout(result,.ctmaTablesIsInit(fit),multi,layout,'study',rowLayout)
  result$layout$sortBy <- sortBy
  result$layout$fontSize <- fontSize
  if(!identical(format,'none')) {
    if(is.null(activeDirectory)) activeDirectory<-if(is.null(models[[1]]$activeDirectory))getwd() else models[[1]]$activeDirectory
    result$files<-.ctmaTablesWrite(result$display,result$notes,format,activeDirectory,file,title,overwrite,layout=layout,fontSize=fontSize)
  }
  result
}
