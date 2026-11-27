# Record the original moderator values alongside the actual TI coding.
.ctmaModeratorCoding <- function(values, predictors, groups, names, centers, scales,
    categories = NULL) {
  lapply(seq_len(ncol(values)), function(i) {
    raw <- values[, i]
    by_study <- split(raw, as.vector(groups))
    study_level <- all(vapply(by_study, function(x) length(unique(x)) == 1L, logical(1)))
    name <- names[i]
    if (length(name) != 1L || is.na(name) || !nzchar(name)) name <- paste("Moderator", i)
    list(name = name, predictors = predictors[[i]], center = centers[[i]],
      scale = scales[[i]], categories = if (is.null(categories)) NULL else categories[[i]],
      reference = if (study_level) raw[!duplicated(groups)] else raw,
      study_level = study_level)
  })
}

.ctmaModeratorProfiles <- function(fit, predictors, moderator, values, profiles, at) {
  if (!is.null(at) && (!is.numeric(at) || any(!is.finite(at)) ||
      is.null(names(at)) || anyDuplicated(names(at)) || !all(names(at) %in% predictors)))
    .ctmaTablesStop("moderatorAt must give finite values named by fitted TI predictors.")
  base <- stats::setNames(numeric(length(predictors)), predictors)
  base[names(at)] <- at
  if (is.null(profiles)) {
    coding <- fit[["moderatorCoding"]]
    if (!length(coding))
      .ctmaTablesStop("This fit has no stored moderator coding. Supply moderatorProfiles with values on the fitted TI scale.")
    index <- if (is.character(moderator)) which(vapply(coding, function(x)
      identical(x$name, moderator), logical(1))) else moderator
    if (!is.numeric(index) || length(index) != 1L || !is.finite(index) ||
        index != floor(index) || index < 1L || index > length(coding))
      .ctmaTablesStop("Select one moderator by its name or position.")
    m <- coding[[index]]
    if (!length(m$predictors) || !all(m$predictors %in% predictors))
      .ctmaTablesStop("The stored moderator coding does not match the fitted TI predictors.")
    if (any(m$predictors %in% names(at)))
      .ctmaTablesStop("Use moderatorValues for the selected moderator and moderatorAt for other predictors.")
    default_values <- is.null(values)
    if (is.null(m$categories)) {
      if (default_values) {
        spread <- stats::sd(m$reference)
        if (!is.finite(spread) || spread <= 0)
          .ctmaTablesStop("The moderator has no estimable SD. Supply moderatorValues explicitly.")
        values <- mean(m$reference) + c(-1, 0, 1)*spread
        labels <- paste0(c("Mean − 1 SD", "Mean", "Mean + 1 SD"),
                         " (", format(signif(values, 7), trim=TRUE), ")")
      } else labels <- as.character(values)
      if (!is.numeric(values) || !length(values) || any(!is.finite(values)))
        .ctmaTablesStop("Continuous moderatorValues must be finite numbers in original units.")
      coded <- matrix((values-m$center)/m$scale, ncol=1L)
    } else {
      if (is.null(values)) values <- m$categories
      if (!length(values) || anyNA(match(values,m$categories)))
        .ctmaTablesStop("moderatorValues contains a category absent from the fitted moderator.")
      labels <- as.character(values)
      coded <- outer(values,m$categories[-1L],`==`) * 1
      coded <- sweep(coded,2L,m$center,`-`)
    }
    profiles <- matrix(rep(base,each=length(values)), nrow=length(values),
                       dimnames=list(labels,predictors))
    profiles[,m$predictors] <- coded
    note <- paste("Moderator", m$name)
    if (is.null(m$categories) && default_values && m$study_level)
      note <- paste0(note, ". Mean and SD give each primary study equal weight.")
    other <- setdiff(predictors,m$predictors)
  } else {
    profiles <- as.matrix(profiles)
    if (!is.numeric(profiles) || !nrow(profiles) || !ncol(profiles) ||
        any(!is.finite(profiles)) || is.null(colnames(profiles)) ||
        anyDuplicated(colnames(profiles)) || !all(colnames(profiles) %in% predictors))
      .ctmaTablesStop("moderatorProfiles must have finite numeric columns named by fitted TI predictors.")
    labels <- rownames(profiles)
    if (is.null(labels)) labels <- paste("Profile",seq_len(nrow(profiles)))
    coded <- matrix(rep(base,each=nrow(profiles)),nrow=nrow(profiles),
                    dimnames=list(labels,predictors))
    other <- setdiff(predictors,colnames(profiles))
    coded[,colnames(profiles)] <- profiles
    profiles <- coded
    note <- "Profiles use the fitted TI predictor scale."
  }
  if (anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels))
    .ctmaTablesStop("Moderator profiles must have distinct labels.")
  list(values=profiles, note=note,
       reference=paste(other,signif(base[other],7),sep=" = ",collapse=", "))
}

.ctmaModeratorStanDrift <- function(fit, raw, predictors) {
  d <- fit$standata
  layout <- .ctmaParMapRawposteriorLayout(d)
  p <- fit$setup$popsetup
  v <- fit$setup$popvalues
  model <- fit$ctstanmodelbase
  registry <- if (exists(".ctMatricesList",asNamespace("ctsem"),inherits=FALSE))
    ".ctMatricesList" else "ctStanMatricesList"
  code <- unname(getFromNamespace(registry,"ctsem")()$all["DRIFT"])
  tform <- getFromNamespace("tform","ctsem")
  stopifnot(ncol(raw) >= max(c(layout$rawpopmeans,layout$tipredeffectparams)))
  out <- array(NA_real_,c(nrow(raw),model$n.latent,model$n.latent))
  fixed <- subset(model$pars, matrix=="DRIFT" & is.na(param) & is.finite(value) &
    row<=model$n.latent & col<=model$n.latent)
  for (i in seq_len(nrow(fixed))) out[,fixed$row[i],fixed$col[i]] <- fixed$value[i]
  for (i in which(p$matrix == code & p$row <= model$n.latent & p$col <= model$n.latent)) {
    source <- if (p$copyrow[i] > 0L) p$copyrow[i] else i
    index <- p$param[source]
    if (index == 0L) {
      out[,p$row[i],p$col[i]] <- v$value[source]
      next
    }
    cell <- raw[,layout$rawpopmeans[index]]
    active <- which(d$TIPREDEFFECTsetup[index,] > 0L)
    if (length(active)) {
      columns <- layout$tipredeffectparams[d$TIPREDEFFECTsetup[index,active]]
      cell <- cell + as.vector(raw[,columns,drop=FALSE] %*% predictors[active])
    }
    out[,p$row[i],p$col[i]] <- tform(cell,p$transform[source],v$multiplier[source],
      v$meanscale[source],v$offset[source],v$inneroffset[source],fit$setup$extratforms)
  }
  out
}

.ctmaTablesModerator <- function(fit, moderator, moderatorValues, moderatorProfiles,
    moderatorAt, moderatorContrast, parameters, statistics, digits, format,
    activeDirectory, file, title, overwrite, tableNote, fontSize, defaultOrder) {
  engine <- .ctmaTablesNative(fit)
  if (is.null(engine)) .ctmaTablesStop("Moderator profiles require a complete pooled fit.")
  if (!is.logical(moderatorContrast) || length(moderatorContrast)!=1L || is.na(moderatorContrast))
    .ctmaTablesStop("moderatorContrast must be TRUE or FALSE.")
  if (!is.numeric(digits) || length(digits)!=1L || !is.finite(digits) ||
      digits!=floor(digits) || digits<0 || digits>12)
    .ctmaTablesStop("digits must be an integer from 0 to 12.")
  if (!is.character(format)||!length(format)||anyNA(format)||
      any(!format %in% c("none","docx","xml","pdf","tex"))||
      ("none" %in% format && length(format)>1L)) .ctmaTablesStop("Unsupported output format.")
  method <- .ctmaPointMethod(fit)
  point_choice <- fit$argumentList$pointEstimate
  if(is.null(point_choice)) point_choice <- "auto"
  raw_point <- NULL
  if (inherits(engine,"ctJuliaFit")) {
    raw_point <- engine$estimate$raw
    raw <- engine$estimate$rawposterior
    if (is.null(raw) || !length(raw)) raw <- matrix(engine$estimate$raw,nrow=1L)
    model <- engine$model_spec$model
    evaluate <- function(x,values=raw) getFromNamespace(".ctBackendPopArrays","ctsem")(
      engine,samples=values,tipreds=x)$pop_DRIFT
  } else {
    if (inherits(engine,"ctmaJuliaFit")) {
      raw_point <- engine$coefficients$raw
      raw <- engine$bayes$draws
      if (is.null(raw)) raw <- engine$transformed_draws
      if (is.null(raw) || !length(raw)) raw <- matrix(engine$coefficients$raw,nrow=1L)
      source <- engine$source_object
      engine <- if (.ctmaJuliaIsUnfittedCoTiMA(source)) .ctmaJuliaGenerateSetup(source) else
        .ctmaJuliaExtractFit(source)
    } else {
      raw_point <- engine$stanfit$rawest
      raw <- engine$stanfit$rawposterior
      if ((is.null(raw) || !length(raw)) && .ctmaTablesStanSampled(engine))
        raw <- getFromNamespace("ctStanRawSamples","ctsem")(engine)
      if (is.null(raw) || !length(raw)) raw <- matrix(engine$stanfit$rawest,nrow=1L)
    }
    model <- engine$ctstanmodelbase
    evaluate <- function(x,values=raw) .ctmaModeratorStanDrift(engine,values,x)
  }
  if (!is.matrix(raw) || !is.numeric(raw) || !length(raw) || any(!is.finite(raw)))
    .ctmaTablesStop("The fit lacks usable joint raw parameter draws or a raw estimate.")
  predictors <- model$TIpredNames
  if (!length(predictors)) .ctmaTablesStop("This model has no fitted TI predictors.")
  profiles <- .ctmaModeratorProfiles(fit,predictors,moderator,moderatorValues,
                                     moderatorProfiles,moderatorAt)
  p <- profiles$values
  if (moderatorContrast && nrow(p)!=2L)
    .ctmaTablesStop("A moderator contrast requires exactly two profiles.")
  d <- model$pars
  d <- d[d$matrix=="DRIFT" & d$row<=model$n.latent & d$col<=model$n.latent,]
  d$parameter <- ifelse(is.na(d$param),paste0("DRIFT[",d$row,",",d$col,"]"),d$param)
  d <- d[.ctmaTablesSelect(d,if(is.null(parameters)) "drift" else parameters),,drop=FALSE]
  if (defaultOrder) d <- d[.ctmaTablesOrder(d),,drop=FALSE]
  stat <- if(is.null(statistics)) "estimate" else tolower(statistics)
  if (length(stat)!=1L || !stat %in% c("mean","median","estimate","mode"))
    .ctmaTablesStop("For moderator profiles, select mean, median, mode, or estimate.")
  point_available <- method%in%c("ML","MAP","optimized") &&
    is.numeric(raw_point) && length(raw_point)==ncol(raw) && all(is.finite(raw_point))
  if(stat=="estimate") stat <- if(point_choice!="mean" && point_available) "mode" else "mean"
  if(stat=="mode" && !point_available) .ctmaTablesStop("Optimized raw values are unavailable for these moderator profiles.")
  points <- if(stat=="mode") lapply(seq_len(nrow(p)),function(i)
    evaluate(p[i,],matrix(as.numeric(raw_point),nrow=1L))) else NULL
  if(length(points)) names(points)<-rownames(p)
  draws <- lapply(seq_len(nrow(p)),function(i) evaluate(p[i,]))
  names(draws) <- rownames(p)
  if (moderatorContrast) {
    contrast <- paste(rownames(p)[2L],"minus",rownames(p)[1L])
    if (contrast %in% names(draws)) .ctmaTablesStop("The contrast label matches a profile label.")
    draws[[contrast]] <- draws[[2L]] - draws[[1L]]
    if(length(points)) points[[contrast]] <- points[[2L]] - points[[1L]]
  }
  display <- data.frame(Parameter=d$parameter)
  data <- lapply(names(draws),function(name) {
    a <- draws[[name]]
    m <- vapply(seq_len(nrow(d)),function(i) a[,d$row[i],d$col[i]],numeric(nrow(raw)))
    m <- matrix(m,nrow=nrow(raw))
    if (any(!is.finite(m))) .ctmaTablesStop("Nonfinite conditional DRIFT values in profile ",name,".")
    estimate <- if(stat=="mode") vapply(seq_len(nrow(d)),function(i)
      points[[name]][1L,d$row[i],d$col[i]],numeric(1)) else
      apply(m,2L,if(stat=="median") stats::median else mean)
    bounds <- if(nrow(m)>1L) apply(m,2L,stats::quantile,probs=c(.025,.975),names=FALSE) else
      matrix(NA_real_,2L,ncol(m))
    f <- function(x) formatC(x,digits=digits,format="f")
    display[[name]] <<- if(nrow(m)>1L)
      paste0(f(estimate)," [",f(bounds[1L,]),", ",f(bounds[2L,]),"]") else f(estimate)
    data.frame(matrix=d$matrix,row=d$row,col=d$col,parameter=d$parameter,profile=name,
               estimate=estimate,lower=bounds[1L,],upper=bounds[2L,])
  })
  notes <- c(profiles$note,if(nzchar(profiles$reference)) paste("Held TI predictors",profiles$reference),
             "DRIFT is reported in the fitted time units.",
             if(nrow(raw)>1L) "Brackets show 95% intervals from joint raw draws." else
               "The saved fit contains a raw point estimate for these profiles.",tableNote)
  notes <- c(notes,if(stat=="mode") paste0("Profile estimates are ",method,
    " point estimates evaluated at the optimized raw parameter vector. Intervals summarize joint uncertainty draws.") else
    paste0("Profile estimates use the ",stat," of transformed joint draws."),.ctmaLegacyPointNotice(fit))
  files <- character()
  if (!identical(format,"none")) {
    if(is.null(activeDirectory)) activeDirectory <- if(is.null(fit$activeDirectory)) getwd() else fit$activeDirectory
    files <- .ctmaTablesWrite(display,notes,format,activeDirectory,file,title,overwrite,
                             layout="wide",fontSize=fontSize)
  }
  structure(list(data=do.call(rbind,data),display=display,notes=notes,files=files,
                 profiles=p,layout=list(layout="wide",fontSize=fontSize)),class="ctmaTables")
}
