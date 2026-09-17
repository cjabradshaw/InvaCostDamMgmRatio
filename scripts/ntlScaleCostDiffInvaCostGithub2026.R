##########################################################################################################
## national-level prediction & extrapolation of costs based on socio-economic traits
## InvaCost Workshop 2022
## January 2024 / updated September 2026
## CJA Bradshaw
##########################################################################################################

library(boot)
library(car)
library(DescTools)
library(dismo)
library(gbm)
library(invacost)
library(lme4)
library(mice)
library(nlme)
library(performance)
library(rcompanion)
library(rworldmap)
library(sf)
library(sjPlot)
library(SpatialEpi)
library(VIM)

# source files
dir.path <- paste(getwd(),"/scripts/",sep="")
source(paste(dir.path,"new_lmer_AIC_tables3.R",sep="")) # change path as required
source(paste(dir.path,"r.squared.R",sep="")) # change path as required
set.seed(101)
dir.create(file.path(getwd(), "out"), recursive=TRUE, showWarnings=FALSE)

# adaptive BRT calibration 
# preserves all existing gbm.step calls below while adapting calibration
# settings to  available sample and retaining diagnostics on  fitted object
`%||%` <- function(x, y) if (is.null(x)) y else x

gbm.step.base <- get("gbm.step", envir = asNamespace("dismo"))
gbm.step.patched <- local({
  gbm.step.fn <- gbm.step.base
  body.lines <- deparse(body(gbm.step.fn), width.cutoff = 500)
  body.lines <- gsub("verbose = FALSE\\)",
                     "verbose = FALSE, n.minobsinnode = n.minobsinnode)",
                     body.lines)
  formals(gbm.step.fn)$n.minobsinnode <- 10L
  body(gbm.step.fn) <- parse(text = paste(body.lines, collapse = "\n"))
  gbm.step.fn
})

gbm.cv.cor <- function(fit) {
  value <- fit$cv.statistics$correlation.mean
  if (length(value) != 1L || !is.finite(value)) -Inf else unname(value)
}

gbm.cv.cor.se <- function(fit) {
  value <- fit$cv.statistics$correlation.se
  if (length(value) != 1L || !is.finite(value)) Inf else unname(value)
}

gbm.cv.jitter <- function(fit) {
  values <- fit$cv.values
  values <- values[is.finite(values)]
  if (length(values) < 7L) return(Inf)
  k <- min(7L, length(values) - 1L + length(values) %% 2L)
  smooth <- stats::runmed(values, k = k, endrule = "median")
  scale.value <- diff(range(values))
  if (!is.finite(scale.value) || scale.value == 0) scale.value <- stats::sd(values)
  if (!is.finite(scale.value) || scale.value == 0) scale.value <- 1
  stats::mad(values - smooth, center = 0, constant = 1) / scale.value
}

gbm.bag.candidates <- function(n, requested = NULL) {
  base <- if (n < 75L) 0.85 else if (n < 150L) 0.75 else if (n < 300L) 0.65 else 0.55
  candidates <- unique(round(c(requested, base, base + 0.1, base - 0.1), 2))
  candidates[is.finite(candidates) & candidates >= 0.5 & candidates <= 0.9]
}

gbm.learning.rate.candidates <- function(requested = NULL) {
  candidates <- unique(c(requested, 0.01, 0.005, 0.001, 0.0005, 0.0002, 0.0001,
                         0.00005, 0.00001))
  candidates[is.finite(candidates) & candidates > 0]
}

gbm.fit.is.better <- function(candidate, candidate.meta, current, current.meta) {
  if (is.null(current)) return(TRUE)
  if (candidate.meta$acceptable != current.meta$acceptable) {
    return(candidate.meta$acceptable)
  }
  candidate.cor <- gbm.cv.cor(candidate)
  current.cor <- gbm.cv.cor(current)
  if (abs(candidate.cor - current.cor) > 0.01) return(candidate.cor > current.cor)
  candidate.jitter <- candidate.meta$jitter
  current.jitter <- current.meta$jitter
  if (abs(candidate.jitter - current.jitter) > 0.02) {
    return(candidate.jitter < current.jitter)
  }
  candidate.se <- gbm.cv.cor.se(candidate)
  current.se <- gbm.cv.cor.se(current)
  if (!isTRUE(all.equal(candidate.se, current.se))) return(candidate.se < current.se)
  abs(candidate$n.trees - candidate.meta$target.min) <
    abs(current$n.trees - current.meta$target.min)
}

adaptive.gbm.step <- function(data, ..., .adaptive.target.min = NULL,
                              .adaptive.refined = FALSE) {
  gbm.args <- list(...)
  gbm.x <- gbm.args$gbm.x
  gbm.y <- gbm.args$gbm.y
  if (is.null(gbm.x) || is.null(gbm.y)) {
    return(do.call(gbm.step.patched, c(list(data = data), gbm.args)))
  }
  n <- nrow(data)
  n.folds <- min(as.integer(gbm.args$n.folds %||% 10L), max(2L, n - 1L))
  max.trees <- gbm.args$max.trees %||% 100000L
  target.min <- .adaptive.target.min %||% if (n < 120L) 500L else 1000L
  tree.complexities <- unique(c(gbm.args$tree.complexity %||% 1L, 1L))
  bag.candidates <- gbm.bag.candidates(n, gbm.args$bag.fraction %||% NULL)
  learning.rates <- gbm.learning.rate.candidates(gbm.args$learning.rate %||% NULL)
  extra.args <- gbm.args[setdiff(names(gbm.args), c(
    "gbm.x", "gbm.y", "weights", "family", "max.trees", "tolerance",
    "learning.rate", "bag.fraction", "tree.complexity", "step.size", "n.folds",
    "silent", "plot.main", "plot.folds", "n.minobsinnode"
  ))]
  best.fit <- NULL
  best.meta <- NULL
  errors <- character()
  attempts <- 0L

  for (tree.complexity in tree.complexities) {
    for (bag.fraction in bag.candidates) {
      step.size <- max(5L, round(n * bag.fraction / 2))
      tolerance <- max(1e-06, 0.001 * sqrt(30 / max(2, n * bag.fraction)))
      n.train <- n * (n.folds - 1L) / n.folds
      n.minobsinnode <- max(2L, min(10L, floor(n.train * bag.fraction / 5)))
      for (learning.rate in learning.rates) {
        attempts <- attempts + 1L
        fit.args <- modifyList(extra.args, list(
          data = data, gbm.x = gbm.x, gbm.y = gbm.y,
          weights = gbm.args$weights %||% NULL,
          family = gbm.args$family %||% "gaussian",
          max.trees = max.trees, tree.complexity = tree.complexity,
          n.folds = n.folds, learning.rate = learning.rate,
          bag.fraction = bag.fraction, step.size = step.size,
          tolerance = tolerance, n.minobsinnode = n.minobsinnode,
          silent = gbm.args$silent %||% FALSE,
          plot.main = gbm.args$plot.main %||% TRUE,
          plot.folds = gbm.args$plot.folds %||% FALSE
        ))
        fit <- tryCatch(
          suppressMessages(suppressWarnings(do.call(gbm.step.patched, fit.args))),
          error = function(error) {
            errors <<- c(errors, conditionMessage(error))
            NULL
          }
        )
        if (is.null(fit)) {
          errors <- c(errors, "gbm.step returned NULL")
          next
        }
        if (length(fit$n.trees) != 1L || !is.finite(fit$n.trees)) {
          errors <- c(errors, "gbm.step returned no finite tree count")
          next
        }
        jitter <- gbm.cv.jitter(fit)
        meta <- list(
          learning.rate = learning.rate, bag.fraction = bag.fraction,
          tree.complexity = tree.complexity, n.folds = n.folds,
          step.size = step.size, tolerance = tolerance,
          n.minobsinnode = n.minobsinnode, target.min = target.min,
          jitter = jitter,
          acceptable = fit$n.trees >= target.min &&
            fit$n.trees < floor(max.trees * 0.95) && jitter <= 0.12
        )
        if (gbm.fit.is.better(fit, meta, best.fit, best.meta)) {
          best.fit <- fit
          best.meta <- meta
        }
        if (meta$acceptable) break
      }
    }
  }
  if (is.null(best.fit)) {
    if (!.adaptive.refined) {
      # Some small samples have no CV improvement under 5-10 folds. Retain the
      # requested conservative learning rate and retry with two folds.
      retry.args <- gbm.args
      retry.args$n.folds <- 2L
      retry.args$tree.complexity <- 1L
      retry.args$bag.fraction <- 0.9
      retry.args$learning.rate <- gbm.args$learning.rate %||% 0.00001
      retry <- do.call(
        adaptive.gbm.step,
        c(list(data = data), retry.args,
          list(.adaptive.target.min = max(200L, round(target.min * 0.5)),
               .adaptive.refined = TRUE))
      )
      retry$adaptive.meta$strategy <- "reduced_folds_simpler_tree"
      return(retry)
    }
    stop("Adaptive gbm.step failed: ", paste(unique(errors), collapse = "; "), call. = FALSE)
  }
  best.fit$adaptive.meta <- c(
    best.meta,
    list(strategy = "default", attempts = attempts, cv.correlation = gbm.cv.cor(best.fit),
         cv.correlation.se = gbm.cv.cor.se(best.fit),
         calibration.errors = unique(errors))
  )
  best.fit
}

gbm.step <- adaptive.gbm.step

AICc <- function(...) {
  models <- list(...)
  vapply(models, function(model) {
    n <- stats::nobs(model)
    k <- attr(stats::logLik(model), "df")
    if (!is.finite(n) || !is.finite(k) || n <= k + 1) return(Inf)
    stats::AIC(model) + (2 * k * (k + 1)) / (n - k - 1)
  }, numeric(1))
}

## functions
in_interval <- function(year, start_year, end_year) {
  !is.na(year) & year >= start_year & year <= end_year
}

delta.IC <- function(x) x - min(x) ## where x is a vector of an IC
weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.IC(AIC.vec); wAIC.vec <- weight.IC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}

# INVACOST
data(invacost)
str(invacost)

#########################
  ## expand costs
  colnames(invacost)
  # eliminating data with no information on starting and ending years
  invacost.db <- invacost[-which(is.na(invacost$Probable_starting_year_adjusted)), ]
  invacost.db <- invacost.db[-which(is.na(invacost.db$Probable_ending_year_adjusted)), ]
  
  # keeping only observed and reliable costs
  invacost.db <- invacost.db[invacost.db$Implementation == "Observed", ]
  invacost.db <- invacost.db[which(invacost.db$Method_reliability == "High"), ]
  
  # eliminating data with no usable cost value
  invacost.db <- invacost.db[-which(is.na(invacost.db$Cost_estimate_per_year_2017_USD_exchange_rate)), ]
  
  invacost.expanded <- expandYearlyCosts(costdb=invacost.db,
                              startcolumn = "Probable_starting_year_adjusted",
                              endcolumn = "Probable_ending_year_adjusted")

  # cycle through countries
  dat.rel <- invacost.expanded
  cntryvec <- sort(unique(trimws(as.character(dat.rel$Official_country))))
  
  cntry.vec <- cntryvec[
    !is.na(cntryvec) &
      nzchar(cntryvec) &
      !grepl("/", cntryvec, fixed = TRUE)
  ]
  
  lcntry <- length(cntry.vec)
  
  stopifnot(
    length(cntry.vec) == length(unique(cntry.vec)),
    all(nzchar(cntry.vec)),
    !anyNA(cntry.vec)
  )
  
  message("number of countries retained: ", lcntry)
  
  excluded.countries <- setdiff(cntryvec, cntry.vec)
  
  write.csv(
    data.frame(country = cntry.vec),
    file.path(getwd(), "out", "countries_retained.csv"),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(country = excluded.countries),
    file.path(getwd(), "out", "countries_excluded.csv"),
    row.names = FALSE
  )
  
  range(invacost.expanded$Impact_year)
  hist(invacost.expanded$Impact_year, xlab="impact year", ylab='frequency', main="")
  length(invacost.expanded$Impact_year[invacost.expanded$Impact_year >= 2000]) / dim(invacost.expanded)[1]
  length(invacost.expanded$Impact_year[invacost.expanded$Impact_year >= 1990]) / dim(invacost.expanded)[1]
  length(invacost.expanded$Impact_year[invacost.expanded$Impact_year >= 1980]) / dim(invacost.expanded)[1]
  
  # discrete three-year bins, with a final two-year bin to retain 2020
  analysis.start.year <- 1980L
  analysis.end.year <- 2020L
  interval.width <- 3L
  
  intervals <- data.frame(
    start_year = seq(
      from = analysis.start.year,
      to = analysis.end.year,
      by = interval.width
    )
  )
  
  intervals$end_year <- pmin(
    intervals$start_year + interval.width - 1L,
    analysis.end.year
  )
  
  st.yr.vec <- intervals$start_year
  en.yr.vec <- intervals$end_year
  
  stopifnot(
    st.yr.vec[1] == analysis.start.year,
    tail(en.yr.vec, 1) == analysis.end.year,
    all(st.yr.vec[-1] == en.yr.vec[-length(en.yr.vec)] + 1L),
    all(en.yr.vec >= st.yr.vec)
  )
  
  print(intervals)
  write.csv(
    intervals,
    file.path(getwd(), "out", "temporal_intervals.csv"),
    row.names = FALSE
  )
  
  iter <- 1000
  itdiv <- iter/10
  
  # country-level point estimates
  cntry.dam.md.boot <- rep(NA_real_, lcntry)
  cntry.mgm.md.boot <- rep(NA_real_, lcntry)
  
  # legacy genus-standardised D:M output
  cntry.dam.mgm.ratio.md.boot <- rep(NA_real_, lcntry)
  
  # raw D:M bootstrap summaries
  cntry.dam.mgm.raw.md.boot <- rep(NA_real_, lcntry)
  cntry.dam.mgm.raw.lo.boot <- rep(NA_real_, lcntry)
  cntry.dam.mgm.raw.up.boot <- rep(NA_real_, lcntry)
  
  # proportion-damage bootstrap summaries
  cntry.pdam.md.boot <- rep(NA_real_, lcntry)
  cntry.pdam.lo.boot <- rep(NA_real_, lcntry)
  cntry.pdam.up.boot <- rep(NA_real_, lcntry)
  
  # proportion-management bootstrap summaries
  cntry.pmgm.md.boot <- rep(NA_real_, lcntry)
  cntry.pmgm.lo.boot <- rep(NA_real_, lcntry)
  cntry.pmgm.up.boot <- rep(NA_real_, lcntry)
  
  # temporal response and data-availability diagnostics
  cntry.dam.mgm.ratio.r.md <- rep(NA_real_, lcntry)
  cntry.usable.intervals <- rep(NA_integer_, lcntry)
  cntry.usable.bootstrap.replicates <- rep(NA_integer_, lcntry)
  
  # genus summaries, retained if still required elsewhere
  cntry.dam.gen.md.boot <- rep(NA_real_, lcntry)
  cntry.mgm.gen.md.boot <- rep(NA_real_, lcntry)
  
  for (c in 1:lcntry) { # 1:lcntry
  
    cntry.dat.rel <- subset(dat.rel, Official_country == cntry.vec[c])
    dim(cntry.dat.rel)
    cntry.dat.rel$Impact_year
  
    # observed costs only
    cntry.obs <- subset(cntry.dat.rel, Implementation == "Observed")
  
    # management costs only
    cntry.obs.mgm <- subset(cntry.obs, Type_of_cost_merged == "Management")
  
    # damage costs only
    cntry.obs.dam <- subset(cntry.obs, Type_of_cost_merged == "Damage")
    
    # Define country × interval × cost-type strata once. Each bootstrap
    # replicate resamples its own observed stratum size, preserving the
    # available data structure without borrowing sample size across types or
    # intervals.
    mgm.by.interval <- lapply(seq_along(st.yr.vec), function(t) {
      cntry.obs.mgm[
        in_interval(cntry.obs.mgm$Impact_year, st.yr.vec[t], en.yr.vec[t]),
        ,
        drop = FALSE
      ]
    })
    
    dam.by.interval <- lapply(seq_along(st.yr.vec), function(t) {
      cntry.obs.dam[
        in_interval(cntry.obs.dam$Impact_year, st.yr.vec[t], en.yr.vec[t]),
        ,
        drop = FALSE
      ]
    })
    
    cntry.mgm.per.md.mat <- cntry.mgm.per.ngen.mat <- cntry.dam.per.md.mat <- cntry.dam.per.ngen.mat <- 
      cntry.mgmPgen.mat <- cntry.damPgen.mat <- matrix(data=NA, ncol=length(st.yr.vec), nrow=iter)
    
    for (t in seq_along(st.yr.vec)) {
      cntry.mgm.per <- mgm.by.interval[[t]]
      cntry.dam.per <- dam.by.interval[[t]]
      cntry.mgml <- nrow(cntry.mgm.per)
      cntry.daml <- nrow(cntry.dam.per)
      
      for (i in seq_len(iter)) {
        if (cntry.mgml > 0L) {
          cntry.mgm.resamp <- cntry.mgm.per[
            sample.int(cntry.mgml, size=cntry.mgml, replace=TRUE),
            ,
            drop = FALSE
          ]
          cntry.mgm.per.md.mat[i,t] <- median(cntry.mgm.resamp$Cost_estimate_per_year_2017_USD_exchange_rate, na.rm=T)
          cntry.mgm.per.ngen.mat[i,t] <- length(table(cntry.mgm.resamp$Genus))
          cntry.mgmPgen.mat[i,t] <- cntry.mgm.per.md.mat[i,t] / cntry.mgm.per.ngen.mat[i,t]
        }
        
        if (cntry.daml > 0L) {
          cntry.dam.resamp <- cntry.dam.per[
            sample.int(cntry.daml, size=cntry.daml, replace=TRUE),
            ,
            drop = FALSE
          ]
          cntry.dam.per.md.mat[i,t] <- median(cntry.dam.resamp$Cost_estimate_per_year_2017_USD_exchange_rate, na.rm=T)
          cntry.dam.per.ngen.mat[i,t] <- length(table(cntry.dam.resamp$Genus))
          cntry.damPgen.mat[i,t] <- cntry.dam.per.md.mat[i,t] / cntry.dam.per.ngen.mat[i,t]
        }
  
        if (i %% itdiv==0) print(i) 
      }
      
      print("###################")
      print(paste(st.yr.vec[t], "to", en.yr.vec[t], sep=" "))
      print("###################")
    }
    
    ###############################################################################
    # COUNTRY-LEVEL BOOTSTRAP SUMMARIES
    #
    # each row is 1 bootstrap replicate.
    # each column is 1 discrete temporal interval.
    #
    # D:M calculated within each bootstrap replicate and interval first.
    # Country-level value for each bootstrap replicate is median
    # across that country's usable temporal intervals
    ###############################################################################
    
    # ---------------------------------------------------------------------------
    # 1. damage-to-management ratio within each bootstrap replicate and interval
    # ---------------------------------------------------------------------------
    
    raw.ratio.boot.mat <- cntry.dam.per.md.mat / cntry.mgm.per.md.mat
    
    # ratios are usable only when both damage and management are finite and
    # management expenditure is strictly positive.
    invalid.ratio <- (
      !is.finite(cntry.dam.per.md.mat) |
        !is.finite(cntry.mgm.per.md.mat) |
        cntry.dam.per.md.mat < 0 |
        cntry.mgm.per.md.mat <= 0 |
        !is.finite(raw.ratio.boot.mat) |
        raw.ratio.boot.mat < 0
    )
    
    raw.ratio.boot.mat[invalid.ratio] <- NA_real_
    
    # ---------------------------------------------------------------------------
    # 2. proportions of total costs attributed to damage and management
    # ---------------------------------------------------------------------------
    
    total.cost.boot.mat <- cntry.dam.per.md.mat + cntry.mgm.per.md.mat
    
    pdam.boot.mat <- cntry.dam.per.md.mat / total.cost.boot.mat
    pmgm.boot.mat <- cntry.mgm.per.md.mat / total.cost.boot.mat
    
    invalid.proportion <- (
      !is.finite(cntry.dam.per.md.mat) |
        !is.finite(cntry.mgm.per.md.mat) |
        cntry.dam.per.md.mat < 0 |
        cntry.mgm.per.md.mat < 0 |
        !is.finite(total.cost.boot.mat) |
        total.cost.boot.mat <= 0
    )
    
    pdam.boot.mat[invalid.proportion] <- NA_real_
    pmgm.boot.mat[invalid.proportion] <- NA_real_
    
    # numerical safeguard against floating-point values slightly outside [0, 1].
    pdam.boot.mat[
      is.finite(pdam.boot.mat) &
        (pdam.boot.mat < 0 | pdam.boot.mat > 1)
    ] <- NA_real_
    
    pmgm.boot.mat[
      is.finite(pmgm.boot.mat) &
        (pmgm.boot.mat < 0 | pmgm.boot.mat > 1)
    ] <- NA_real_
    
    # ---------------------------------------------------------------------------
    # 3. safe row-median helper
    #
    # apply(..., median, na.rm = TRUE) returns Inf or warnings for rows containing
    # no finite values. This helper explicitly returns NA for such rows.
    # ---------------------------------------------------------------------------
    
    row_median_finite <- function(x) {
      x <- x[is.finite(x)]
      
      if (length(x) == 0L) {
        return(NA_real_)
      }
      
      stats::median(x)
    }
    
    country.dam.boot <- apply(
      cntry.dam.per.md.mat,
      MARGIN = 1L,
      FUN = row_median_finite
    )
    country.mgm.boot <- apply(
      cntry.mgm.per.md.mat,
      MARGIN = 1L,
      FUN = row_median_finite
    )
    country.dam.gen.boot <- apply(
      cntry.dam.per.ngen.mat,
      MARGIN = 1L,
      FUN = row_median_finite
    )
    country.mgm.gen.boot <- apply(
      cntry.mgm.per.ngen.mat,
      MARGIN = 1L,
      FUN = row_median_finite
    )
    
    genus.ratio.boot.mat <- cntry.damPgen.mat / cntry.mgmPgen.mat
    genus.ratio.boot.mat[
      !is.finite(genus.ratio.boot.mat) | cntry.mgmPgen.mat <= 0
    ] <- NA_real_
    country.genus.ratio.boot <- apply(
      genus.ratio.boot.mat,
      MARGIN = 1L,
      FUN = row_median_finite
    )
    
    # Rate is calculated only between temporally adjacent usable intervals.
    # Taking differences after omitting missing intervals would incorrectly
    # treat non-adjacent observations as consecutive periods.
    ratio.rate.boot.mat <- matrix(
      NA_real_,
      nrow = iter,
      ncol = max(0L, ncol(raw.ratio.boot.mat) - 1L)
    )
    if (ncol(raw.ratio.boot.mat) > 1L) {
      previous.ratio <- raw.ratio.boot.mat[, -ncol(raw.ratio.boot.mat), drop=FALSE]
      next.ratio <- raw.ratio.boot.mat[, -1L, drop=FALSE]
      valid.rate <- is.finite(previous.ratio) & is.finite(next.ratio) &
        previous.ratio > 0 & next.ratio > 0
      ratio.rate.boot.mat[valid.rate] <- log10(
        next.ratio[valid.rate] / previous.ratio[valid.rate]
      )
    }
    country.ratio.rate.boot <- apply(
      ratio.rate.boot.mat,
      MARGIN = 1L,
      FUN = row_median_finite
    )
    
    # ---------------------------------------------------------------------------
    # 4. one country-level value per bootstrap replicate
    #
    # The median is calculated across usable temporal intervals separately for
    # every bootstrap replicate.
    # ---------------------------------------------------------------------------
    
    country.ratio.boot <- apply(
      raw.ratio.boot.mat,
      MARGIN = 1L,
      FUN = row_median_finite
    )
    
    country.pdam.boot <- apply(
      pdam.boot.mat,
      MARGIN = 1L,
      FUN = row_median_finite
    )
    
    country.pmgm.boot <- apply(
      pmgm.boot.mat,
      MARGIN = 1L,
      FUN = row_median_finite
    )

    # Auxiliary country summaries required by the later exported tables.
    if (any(is.finite(country.dam.boot))) {
      cntry.dam.md.boot[c] <- stats::median(country.dam.boot, na.rm=TRUE)
    }
    if (any(is.finite(country.mgm.boot))) {
      cntry.mgm.md.boot[c] <- stats::median(country.mgm.boot, na.rm=TRUE)
    }
    if (any(is.finite(country.dam.gen.boot))) {
      cntry.dam.gen.md.boot[c] <- stats::median(country.dam.gen.boot, na.rm=TRUE)
    }
    if (any(is.finite(country.mgm.gen.boot))) {
      cntry.mgm.gen.md.boot[c] <- stats::median(country.mgm.gen.boot, na.rm=TRUE)
    }
    if (any(is.finite(country.genus.ratio.boot))) {
      cntry.dam.mgm.ratio.md.boot[c] <- stats::median(
        country.genus.ratio.boot,
        na.rm=TRUE
      )
    }
    if (any(is.finite(country.ratio.rate.boot))) {
      cntry.dam.mgm.ratio.r.md[c] <- stats::median(
        country.ratio.rate.boot,
        na.rm=TRUE
      )
    }
    
    # Retain only finite bootstrap estimates.
    country.ratio.boot[!is.finite(country.ratio.boot)] <- NA_real_
    country.pdam.boot[!is.finite(country.pdam.boot)] <- NA_real_
    country.pmgm.boot[!is.finite(country.pmgm.boot)] <- NA_real_
    
    # ---------------------------------------------------------------------------
    # 5. record data-availability diagnostics
    # ---------------------------------------------------------------------------
    
    # Number of intervals for which at least one bootstrap replicate produced
    # a valid D:M estimate.
    cntry.usable.intervals[c] <- sum(
      colSums(is.finite(raw.ratio.boot.mat)) > 0L
    )
    
    # Number of bootstrap replicates producing a valid country-level D:M estimate.
    cntry.usable.bootstrap.replicates[c] <- sum(
      is.finite(country.ratio.boot)
    )
    
    # ---------------------------------------------------------------------------
    # 6. summarise the country-level D:M bootstrap distribution
    # ---------------------------------------------------------------------------
    
    if (sum(is.finite(country.ratio.boot)) > 0L) {
      
      cntry.dam.mgm.raw.md.boot[c] <- stats::median(
        country.ratio.boot,
        na.rm = TRUE
      )
      
      cntry.dam.mgm.raw.lo.boot[c] <- stats::quantile(
        country.ratio.boot,
        probs = 0.025,
        na.rm = TRUE,
        names = FALSE,
        type = 7
      )
      
      cntry.dam.mgm.raw.up.boot[c] <- stats::quantile(
        country.ratio.boot,
        probs = 0.975,
        na.rm = TRUE,
        names = FALSE,
        type = 7
      )
    }
    
    # ---------------------------------------------------------------------------
    # 7. summarise the country-level pdam bootstrap distribution
    # ---------------------------------------------------------------------------
    
    if (sum(is.finite(country.pdam.boot)) > 0L) {
      
      cntry.pdam.md.boot[c] <- stats::median(
        country.pdam.boot,
        na.rm = TRUE
      )
      
      cntry.pdam.lo.boot[c] <- stats::quantile(
        country.pdam.boot,
        probs = 0.025,
        na.rm = TRUE,
        names = FALSE,
        type = 7
      )
      
      cntry.pdam.up.boot[c] <- stats::quantile(
        country.pdam.boot,
        probs = 0.975,
        na.rm = TRUE,
        names = FALSE,
        type = 7
      )
    }
    
    # ---------------------------------------------------------------------------
    # 8. summarise the country-level pmgm bootstrap distribution
    # ---------------------------------------------------------------------------
    
    if (sum(is.finite(country.pmgm.boot)) > 0L) {
      
      cntry.pmgm.md.boot[c] <- stats::median(
        country.pmgm.boot,
        na.rm = TRUE
      )
      
      cntry.pmgm.lo.boot[c] <- stats::quantile(
        country.pmgm.boot,
        probs = 0.025,
        na.rm = TRUE,
        names = FALSE,
        type = 7
      )
      
      cntry.pmgm.up.boot[c] <- stats::quantile(
        country.pmgm.boot,
        probs = 0.975,
        na.rm = TRUE,
        names = FALSE,
        type = 7
      )
    }
    
    # ---------------------------------------------------------------------------
    # 9. internal consistency checks
    # ---------------------------------------------------------------------------
    
    if (
      is.finite(cntry.pdam.md.boot[c]) &&
      is.finite(cntry.pmgm.md.boot[c])
    ) {
      stopifnot(
        abs(
          cntry.pdam.md.boot[c] +
            cntry.pmgm.md.boot[c] -
            1
        ) < 1e-8
      )
    }
    
    if (is.finite(cntry.dam.mgm.raw.md.boot[c])) {
      stopifnot(
        cntry.dam.mgm.raw.md.boot[c] >= 0,
        cntry.dam.mgm.raw.lo.boot[c] <=
          cntry.dam.mgm.raw.md.boot[c],
        cntry.dam.mgm.raw.up.boot[c] >=
          cntry.dam.mgm.raw.md.boot[c]
      )
    }  
    
    print("###################")
    print(paste(cntry.vec[c], " (", round(c/lcntry*100, 0), "% complete)", sep=""))
    print("###################")
  } # end c
  
ratio.out <- data.frame(
  country = cntry.vec,
  
  # legacy genus-standardised value, if retained
  damPgen.mgmPgen = cntry.dam.mgm.ratio.md.boot,
  
  # principal raw D:M estimate and bootstrap interval
  dam.mgm = cntry.dam.mgm.raw.md.boot,
  dam.mgm.lo = cntry.dam.mgm.raw.lo.boot,
  dam.mgm.up = cntry.dam.mgm.raw.up.boot,
  
  # temporal response
  ratio.r = cntry.dam.mgm.ratio.r.md,
  
  # proportion damage and bootstrap interval
  pdam = cntry.pdam.md.boot,
  pdam.lo = cntry.pdam.lo.boot,
  pdam.up = cntry.pdam.up.boot,
  
  # proportion management and bootstrap interval
  pmgm = cntry.pmgm.md.boot,
  pmgm.lo = cntry.pmgm.lo.boot,
  pmgm.up = cntry.pmgm.up.boot,
  
  # data-availability diagnostics
  usable.intervals = cntry.usable.intervals,
  usable.bootstrap.replicates =
    cntry.usable.bootstrap.replicates,
  
  stringsAsFactors = FALSE
)

## validation
stopifnot(
  nrow(ratio.out) == length(cntry.vec),
  !anyDuplicated(ratio.out$country)
)

valid.proportions <- with(
  ratio.out,
  is.finite(pdam) & is.finite(pmgm)
)

if (any(valid.proportions)) {
  stopifnot(
    all(
      abs(
        ratio.out$pdam[valid.proportions] +
          ratio.out$pmgm[valid.proportions] -
          1
      ) < 1e-8
    )
  )
}

dir.path <- paste(getwd(),"/out/",sep="")
write.csv(ratio.out, paste(dir.path,"ratioOut.csv",sep=""))
hist(log10(ratio.out$dam.mgm),main="",xlab="log10 damage:management")
hist((ratio.out$ratio.r),main="",xlab="damage:management r")
hist(logit(ratio.out$pdam),main="",xlab="damage:management r")
hist(logit(ratio.out$pmgm),main="",xlab="damage:management r")

# mgm by dam relationship (across countries)
dam.mgm.out <- data.frame(cntry.vec, cntry.dam.md.boot, cntry.mgm.md.boot)
colnames(dam.mgm.out) <- c("country","dam","mgm")
plot(log10(dam.mgm.out$dam), log10(dam.mgm.out$mgm), pch=19)
fit.dam.mgm <- lm(log10(mgm) ~ log10(dam), data=dam.mgm.out)
summary(fit.dam.mgm)
abline(fit.dam.mgm, lty=2, col="red")
linreg.ER(log10(dam.mgm.out$mgm), log10(dam.mgm.out$dam))

check_model(fit.dam.mgm)
plot_model(fit.dam.mgm, show.values=T, vline.color = "purple")

# dam by number of genera
dam.ngen.out <- data.frame(cntry.vec, cntry.dam.gen.md.boot, cntry.dam.md.boot)
colnames(dam.ngen.out) <- c("country", "ngen", "dam")
dam.nggen.out.noNA <- na.omit(dam.ngen.out)
head(dam.nggen.out.noNA)
plot(log(dam.nggen.out.noNA$ngen), (dam.nggen.out.noNA$dam), pch=19, xlab="number of genera", ylab="log annual damage cost")
fit <- lm(log(dam.nggen.out.noNA$dam) ~ (dam.nggen.out.noNA$ngen))
summary(fit)
abline(fit, lty=2, col="red")

# mgm by number of genera
mgm.ngen.out <- data.frame(cntry.vec, cntry.mgm.gen.md.boot, cntry.mgm.md.boot)
colnames(mgm.ngen.out) <- c("country", "ngen", "mgm")
head(mgm.ngen.out)
plot(log10(mgm.ngen.out$ngen), log10(mgm.ngen.out$mgm), pch=19, xlab="log10 number of genera", ylab="log10 annual management cost")
fit <- lm(log10(mgm.ngen.out$mgm) ~ log10(mgm.ngen.out$ngen))
summary(fit)
abline(fit, lty=2, col="red")

# import cntry.code
dir.path <- paste(getwd(),"/data/",sep="")
cont.cntry <- read.csv(paste(dir.path,"continent.countryINVACOST.csv",sep=""), header=T)

# merge with cntry.code
ratio.cntry <- merge(ratio.out, cont.cntry, by="country")
head(ratio.cntry)

dam.mgm.cntry.out <- merge(dam.mgm.out, cont.cntry, by="country")
table(dam.mgm.cntry.out$cont)
dam.mgm.cntry.out$reg2 <- ifelse(dam.mgm.cntry.out$cont == "NAM" | dam.mgm.cntry.out$cont == "CAR", "NAMCAR", dam.mgm.cntry.out$cont)
dam.mgm.cntry.out$reg2 <- ifelse(dam.mgm.cntry.out$cont == "EUR" | dam.mgm.cntry.out$cont == "ME", "EURME", dam.mgm.cntry.out$reg2)
dam.mgm.cntry.out$reg2 <- ifelse(dam.mgm.cntry.out$cont == "ASIA" | dam.mgm.cntry.out$cont == "OC", "ASIAOC", dam.mgm.cntry.out$reg2)
table(dam.mgm.cntry.out$reg2)

dir.path <- paste(getwd(),"/out/",sep="")
write.csv(dam.mgm.cntry.out, paste(dir.path,"damMgmOut.csv",sep=""))

dam.ngen.cntry.out <- merge(dam.ngen.out, cont.cntry, by="country")
table(dam.ngen.cntry.out$cont)
dam.ngen.cntry.out$reg2 <- ifelse(dam.ngen.cntry.out$cont == "NAM" | dam.ngen.cntry.out$cont == "CAR", "NAMCAR", dam.ngen.cntry.out$cont)
dam.ngen.cntry.out$reg2 <- ifelse(dam.ngen.cntry.out$cont == "EUR" | dam.ngen.cntry.out$cont == "ME", "EURME", dam.ngen.cntry.out$reg2)
dam.ngen.cntry.out$reg2 <- ifelse(dam.ngen.cntry.out$cont == "ASIA" | dam.ngen.cntry.out$cont == "OC", "ASIAOC", dam.ngen.cntry.out$reg2)
table(dam.ngen.cntry.out$reg2)

write.csv(dam.ngen.cntry.out, paste(dir.path,"damGenOut.csv", sep=""))

mgm.ngen.cntry.out <- merge(mgm.ngen.out, cont.cntry, by="country")
table(mgm.ngen.cntry.out$cont)
mgm.ngen.cntry.out$reg2 <- ifelse(mgm.ngen.cntry.out$cont == "NAM" | mgm.ngen.cntry.out$cont == "CAR", "NAMCAR", mgm.ngen.cntry.out$cont)
mgm.ngen.cntry.out$reg2 <- ifelse(mgm.ngen.cntry.out$cont == "EUR" | mgm.ngen.cntry.out$cont == "ME", "EURME", mgm.ngen.cntry.out$reg2)
mgm.ngen.cntry.out$reg2 <- ifelse(mgm.ngen.cntry.out$cont == "ASIA" | mgm.ngen.cntry.out$cont == "OC", "ASIAOC", mgm.ngen.cntry.out$reg2)
table(mgm.ngen.cntry.out$reg2)

dir.path <- paste(getwd(),"/out/",sep="")
write.csv(mgm.ngen.cntry.out, paste(dir.path,"mgmGenOut.csv", sep=""))


# GDP per capita
dir.path <- paste(getwd(),"/data/",sep="")
gdp <- read.csv(paste(dir.path,"GDPpc.csv", sep=""))
head(gdp)
ratio.gdp <- merge(ratio.cntry, gdp, by="cntry.code", all = T)
gdp.values <- trimws(as.character(ratio.gdp$GDPpcRecent))
gdp.missing <- is.na(gdp.values) | gdp.values %in% c("", "#DIV/0!")
gdp.numeric <- suppressWarnings(as.numeric(gdp.values))
if (any(!gdp.missing & is.na(gdp.numeric))) {
  stop("GDPpcRecent contains unrecognised non-numeric values.", call.=FALSE)
}
ratio.gdp$GDPpcRecent <- gdp.numeric
stopifnot(is.numeric(ratio.gdp$GDPpcRecent))
head(ratio.gdp)

plot(log10(as.numeric(ratio.gdp$GDPpcRecent)), log10(ratio.gdp$dam.mgm), pch=19, ylab="log10 damage:management", xlab="log10 GDP per capita")
fit <- lm(log10(ratio.gdp$dam.mgm) ~ log10(ratio.gdp$GDPpcRecent))
summary(fit)
abline(fit,lty=2, col="red")

plot(log10(ratio.gdp$GDPpcRecent), (ratio.gdp$ratio.r), pch=19, ylab="damage:management r", xlab="log10 GDP per capita")
fit <- lm((ratio.gdp$ratio.r) ~ log10(ratio.gdp$GDPpcRecent))
summary(fit)
abline(fit,lty=2, col="red")

# by continent
table(ratio.gdp$cont)
EUR.dat <- subset(ratio.gdp, cont=="EUR")
plot(log10(EUR.dat$GDPpcRecent), log10(EUR.dat$dam.mgm), pch=19, ylab="log10 damage:management", xlab="log10 GDP per capita")
fit <- lm(log10(EUR.dat$dam.mgm) ~ log10(EUR.dat$GDPpcRecent))
abline(fit,lty=2, col="red")

ASIA.dat <- subset(ratio.gdp, cont=="ASIA")
plot(log10(ASIA.dat$GDPpcRecent), log10(ASIA.dat$dam.mgm), pch=19, ylab="log10 damage:management", xlab="log10 GDP per capita")
fit <- lm(log10(ASIA.dat$dam.mgm) ~ log10(ASIA.dat$GDPpcRecent))
abline(fit,lty=2, col="red")

AFR.dat <- subset(ratio.gdp, cont=="AFR")
plot(log10(AFR.dat$GDPpcRecent), log10(AFR.dat$dam.mgm), pch=19, ylab="log10 damage:management", xlab="log10 GDP per capita")
fit <- lm(log10(AFR.dat$dam.mgm) ~ log10(AFR.dat$GDPpcRecent))
abline(fit,lty=2, col="red")

NAM.dat <- subset(ratio.gdp, cont=="NAM")
plot(log10(NAM.dat$GDPpcRecent), log10(NAM.dat$dam.mgm), pch=19, ylab="log10 damage:management", xlab="log10 GDP per capita")
fit <- lm(log10(NAM.dat$dam.mgm) ~ log10(NAM.dat$GDPpcRecent))
abline(fit,lty=2, col="red")

SA.dat <- subset(ratio.gdp, cont=="SA")
plot(log10(SA.dat$GDPpcRecent), log10(SA.dat$dam.mgm), pch=19, ylab="log10 damage:management", xlab="log10 GDP per capita")
fit <- lm(log10(SA.dat$dam.mgm) ~ log10(SA.dat$GDPpcRecent))
abline(fit,lty=2, col="red")

OC.dat <- subset(ratio.gdp, cont=="OC")
plot(log10(OC.dat$GDPpcRecent), log10(OC.dat$dam.mgm), pch=19, ylab="log10 damage:management", xlab="log10 GDP per capita")
fit <- lm(log10(OC.dat$dam.mgm) ~ log10(OC.dat$GDPpcRecent))
abline(fit,lty=2, col="red")

## FAO proportion of GDP devoted to agriculture
## Agriculture, forestry, and fishing, value added (% of GDP)
## https://data.worldbank.org/indicator/NV.AGR.TOTL.ZS
# Agriculture, forestry, and fishing corresponds to ISIC divisions 1-3 and includes forestry, hunting, and fishing, as well as
# cultivation of crops and livestock production. Value added is the net output of a sector after adding up all outputs and subtracting
# intermediate inputs. It is calculated without making deductions for depreciation of fabricated assets or depletion and degradation
# of natural resources. The origin of value added is determined by the International Standard Industrial Classification (ISIC),
# revision 4. Note: For VAB countries, gross value added at factor cost is used as the denominator.
FAOag <- read.csv(paste(dir.path,"faoag.csv",sep=""))
VAag.all <- subset(FAOag, Item == 'Value Added (Agriculture, Forestry and Fishing)' & ElementCode == 6139)
head(VAag.all)

# cycle through countries for most recent value
VAag.curr <- VAag.all[1,]
cntry.vec <- names(table(VAag.all$Area))
for (c in 1:length(cntry.vec)) {
  cntry.it <- subset(VAag.all, Area == cntry.vec[c])
  VAag.curr <- rbind(VAag.curr, cntry.it[which(cntry.it$Year == max(cntry.it$Year, na.rm=T)),])
}
VAag.curr <- VAag.curr[-1,]
head(VAag.curr)
VAag1 <- VAag.curr[,c(4,10,12)]
colnames(VAag1) <- c("country", "year", "VAag")
head(VAag1)

fao.cntry.code <- read.csv(paste(dir.path,"fao.cntry.code.csv", sep=""))
head(fao.cntry.code)
VAag <- merge(VAag1, fao.cntry.code, by="country", all=T)
head(VAag)

# Transparency International
# Brenton-Rule et al. 2016 - corruption (no relationship with # invasive species)
# https://royalsocietypublishing.org/doi/full/10.1098/rspb.2016.0901
# Corruption Perception index https://www.transparency.org/en/cpi/2021
cpi <- read.csv(paste(dir.path,"CPI.csv",sep=""))
head(cpi)

# gov expenditure on education (% GDP)
# https://data.worldbank.org/indicator/SE.XPD.TOTL.GD.ZS
expedu <- read.csv(paste(dir.path,"govexpedu.csv",sep=""))
head(expedu)

# Global Health Security Index
# https://www.ghsindex.org/report-model/
ghsi <- read.csv(paste(dir.path,"GHSI2022.csv",sep=""))
head(ghsi)

# % agricultural land
# https://data.worldbank.org/indicator/AG.LND.AGRI.ZS
agrlnd <- read.csv(paste(dir.path,"pcAgrLand.csv",sep=""))
head(agrlnd)

# imports of goods and services per capita (mean of last 5 years)
# https://data.worldbank.org/indicator/NE.IMP.GNFS.CD
igs <- read.csv(paste(dir.path,"importGS.csv",sep=""))
head(igs)
pop21 <- read.csv(paste(dir.path,"pop2021.csv", sep="")) # https://data.worldbank.org/indicator/SP.POP.TOTL
head(pop21)
igspc <- merge(igs, pop21, by="cntry.code", all=T)
head(igspc)
igspc$igspc <- igspc$igs/igspc$popN
head(igspc)

# scientific & technical journal articles (per capita)
# https://data.worldbank.org/indicator/IP.JRN.ARTC.SC
stja <- read.csv(paste(dir.path,"stjarticles.csv",sep=""))
head(stja)
stjapc <- merge(stja, pop21, by="cntry.code", all=T)
head(stjapc)
stjapc$stjapc <- stjapc$stja/stjapc$popN
head(stjapc)

ratio.gdpcpi <- merge(ratio.gdp, cpi, by="cntry.code", all=T)
head(ratio.gdpcpi)

ratio.full1 <- merge(ratio.gdpcpi, expedu, by="cntry.code", all=T)
head(ratio.full1)

ratio.full2 <- merge(ratio.full1, VAag, by="cntry.code", all=T)
head(ratio.full2)

ratio.full3 <- merge(ratio.full2, ghsi, by="cntry.code", all=T)
head(ratio.full3)

ratio.full4 <- merge(ratio.full3, agrlnd, by="cntry.code", all=T)
head(ratio.full4)

ratio.full5 <- merge(ratio.full4, igspc, by="cntry.code", all=T)
head(ratio.full5)

ratio.full6 <- merge(ratio.full5, stjapc, by="cntry.code", all=T)
head(ratio.full6)

ratio.full <- ratio.full6
head(ratio.full)

dim(ratio.full)
head(ratio.full)
length(table(ratio.full$cntry.code))
ratio.full <- ratio.full[is.na(ratio.full$cntry.code)==F,]
dim(ratio.full)
head(ratio.full)
length(table(ratio.full$cntry.code))


plot((ratio.full$CPI), log10(ratio.full$dam.mgm), pch=19, ylab="log10 damage:management", xlab="corruption index")
fit <- lm(log10(ratio.full$dam.mgm) ~ (ratio.full$CPI))
abline(fit,lty=2, col="red")

plot((ratio.full$govexpedu), log10(ratio.full$dam.mgm), pch=19, ylab="log10 damage:management", xlab="gov education expend %GDP")
fit <- lm(log10(ratio.full$dam.mgm) ~ (ratio.full$govexpedu))
abline(fit,lty=2, col="red")

plot((ratio.full$VAag), log10(ratio.full$dam.mgm), pch=19, ylab="log10 damage:management", xlab="value added agriculture/fisheries/forestry %GDP")
fit <- lm(log10(ratio.full$dam.mgm) ~ (ratio.full$VAag))
abline(fit,lty=2, col="red")

plot((ratio.full$GHSI), log10(ratio.full$dam.mgm), pch=19, ylab="log10 damage:management", xlab="global health security index")
fit <- lm(log10(ratio.full$dam.mgm) ~ (ratio.full$GHSI))
abline(fit,lty=2, col="red")

plot((ratio.full$pcAgrLand), log10(ratio.full$dam.mgm), pch=19, ylab="log10 damage:management", xlab="% agricultural land")
fit <- lm(log10(ratio.full$dam.mgm) ~ (ratio.full$pcAgrLand))
abline(fit,lty=2, col="red")

plot(log10(ratio.full$igspc), log10(ratio.full$dam.mgm), pch=19, ylab="log10 damage:management", xlab="log10 pc import goods & services")
fit <- lm(log10(ratio.full$dam.mgm) ~ log10(ratio.full$igspc))
abline(fit,lty=2, col="red")

plot(log10(ratio.full$stjapc), log10(ratio.full$dam.mgm), pch=19, ylab="log10 damage:management", xlab="log10 pc sci/tech articles")
fit <- lm(log10(ratio.full$dam.mgm) ~ log10(ratio.full$stjapc))
abline(fit,lty=2, col="red")

plot(log10(ratio.full$GDPpcRecent), logit(ratio.full$pdam), pch=19, ylab="logit proportion damage", xlab="log10 GDP per capita")
fit <- lm(logit(ratio.full$pdam) ~ log10(ratio.full$GDPpcRecent))
summary(fit)
abline(fit,lty=2, col="red")
plot(log10(ratio.full$GDPpcRecent), logit(ratio.full$pmgm), pch=19, ylab="logit proportion management", xlab="log10 GDP per capita")
fit <- lm(logit(ratio.full$pmgm) ~ log10(ratio.full$GDPpcRecent))
abline(fit,lty=2, col="red")

plot((ratio.full$CPI), logit(ratio.full$pdam), pch=19, ylab="logit proportion damage", xlab="corruption index")
fit <- lm(logit(ratio.full$pdam) ~ (ratio.full$CPI))
abline(fit,lty=2, col="red")
plot((ratio.full$CPI), logit(ratio.full$pmgm), pch=19, ylab="logit proportion management", xlab="corruption index")
fit <- lm(logit(ratio.full$pmgm) ~ (ratio.full$CPI))
abline(fit,lty=2, col="red")

plot((ratio.full$govexpedu), logit(ratio.full$pdam), pch=19, ylab="logit proportion damage", xlab="gov education expend %GDP")
fit <- lm(logit(ratio.full$pdam) ~ (ratio.full$govexpedu))
abline(fit,lty=2, col="red")
plot((ratio.full$govexpedu), logit(ratio.full$pmgm), pch=19, ylab="logit proportion management", xlab="gov education expend %GDP")
fit <- lm(logit(ratio.full$pmgm) ~ (ratio.full$govexpedu))
abline(fit,lty=2, col="red")

plot((ratio.full$VAag/100), logit(ratio.full$pdam), pch=19, ylab="logit proportion damage", xlab="value added primary production %GDP")
fit <- lm(logit(ratio.full$pdam) ~ (ratio.full$VAag))
abline(fit,lty=2, col="red")
plot((ratio.full$VAag), logit(ratio.full$pmgm), pch=19, ylab="logit proportion management", xlab="value added primary production %GDP")
fit <- lm(logit(ratio.full$pmgm) ~ (ratio.full$VAag))
abline(fit,lty=2, col="red")

plot((ratio.full$GHSI), logit(ratio.full$pdam), pch=19, ylab="logit proportion damage", xlab="global health security index")
fit <- lm(logit(ratio.full$pdam) ~ (ratio.full$GHSI))
abline(fit,lty=2, col="red")
plot((ratio.full$GHSI), logit(ratio.full$pmgm), pch=19, ylab="logit proportion management", xlab="global health security index")
fit <- lm(logit(ratio.full$pmgm) ~ (ratio.full$GHSI))
abline(fit,lty=2, col="red")

plot(logit(ratio.full$pcAgrLand/100), logit(ratio.full$pdam), pch=19, ylab="logit proportion damage", xlab="logit % agricultural land")
fit <- lm(logit(ratio.full$pdam) ~ logit(ratio.full$pcAgrLand/100))
abline(fit,lty=2, col="red")
plot(logit(ratio.full$pcAgrLand/100), logit(ratio.full$pmgm), pch=19, ylab="logit proportion management", xlab="logit % agricultural land")
fit <- lm(logit(ratio.full$pmgm) ~ logit(ratio.full$pcAgrLand/100))
abline(fit,lty=2, col="red")

plot(log10(ratio.full$igspc), logit(ratio.full$pdam), pch=19, ylab="logit proportion damage", xlab="log10 pc import goods & services")
fit <- lm(logit(ratio.full$pdam) ~ log10(ratio.full$igspc))
summary(fit)
abline(fit,lty=2, col="red")
plot(log10(ratio.full$igspc), logit(ratio.full$pmgm), pch=19, ylab="logit proportion management", xlab="log10 pc import goods & services")
fit <- lm(logit(ratio.full$pmgm) ~ log10(ratio.full$igspc))
abline(fit,lty=2, col="red")

plot(log10(ratio.full$stjapc), logit(ratio.full$pdam), pch=19, ylab="logit proportion damage", xlab="log10 pc sci/tech articles")
fit <- lm(logit(ratio.full$pdam) ~ log10(ratio.full$stjapc))
summary(fit)
abline(fit,lty=2, col="red")
plot(log10(ratio.full$stjapc), logit(ratio.full$pmgm), pch=19, ylab="logit proportion management", xlab="log10 pc sci/tech articles")
fit <- lm(logit(ratio.full$pmgm) ~ log10(ratio.full$stjapc))
abline(fit,lty=2, col="red")


final.dat.raw <- data.frame(ratio.full$GDPpcRecent, ratio.full$CPI, ratio.full$GHSI, ratio.full$igspc, ratio.full$pcAgrLand,
                            ratio.full$VAag, ratio.full$govexpedu, ratio.full$stjapc)
colnames(final.dat.raw) <- c("gdp", "cpi", "ghsi", "igs", "agrL", "VAag", "govexpedu", "stja")
head(final.dat.raw)
tail(final.dat.raw)
dim(final.dat.raw)
dim(na.omit(final.dat.raw))

aggr_plot <- aggr(final.dat.raw, col=c('navyblue','red'), numbers=FALSE, sortVars=TRUE, labels=names(final.dat.raw), cex.axis=.7, gap=3, ylab=c("histogram of missing data","pattern"))
marginplot(final.dat.raw[c(1,2)])
aggr_plot$missings
missings <- data.frame(aggr_plot$missings$Variable, round(100*aggr_plot$missings$Count/dim(final.dat.raw)[1], 1))
colnames(missings) <- c("var", "pcMiss")
missings

final.dat.rawRatio <- data.frame(ratio.full$dam.mgm, final.dat.raw)
colnames(final.dat.rawRatio)[1] <- "dam.mgm"
head(final.dat.rawRatio)
final.dat.rawRatio <- final.dat.rawRatio[is.na(final.dat.rawRatio$dam.mgm)==F,]
head(final.dat.rawRatio)
dim(final.dat.rawRatio)
aggr_plot2 <- aggr(final.dat.rawRatio, col=c('navyblue','red'), numbers=FALSE, sortVars=TRUE, labels=names(final.dat.raw), cex.axis=.7, gap=3, ylab=c("histogram of missing data","pattern"))


########################
# multiple imputation
########################
imputation.keep <- rowSums(!is.na(final.dat.raw)) > 0
imputation.data <- cbind(
  dam.mgm=ratio.full$dam.mgm, pdam=ratio.full$pdam, ratio.r=ratio.full$ratio.r,
  final.dat.raw
)
imputation.data <- imputation.data[imputation.keep, , drop=FALSE]
if (!all(vapply(final.dat.raw, is.numeric, logical(1)))) {
  stop("All imputed predictors must be numeric.", call.=FALSE)
}
if (any(vapply(final.dat.raw, function(x) any(!is.na(x) & !is.finite(x)), logical(1)))) {
  stop("Imputed predictors contain non-finite values.", call.=FALSE)
}
imputation.method <- make.method(imputation.data)
imputation.method[c("dam.mgm", "pdam", "ratio.r")] <- ""
imputation.predictors <- make.predictorMatrix(imputation.data)
imputation.predictors[c("dam.mgm", "pdam", "ratio.r"), ] <- 0

# mice consistently identifies this predictor pair as collinear when imputing
# value-added agriculture, so exclude it
imputation.predictors["VAag", "stja"] <- 0
md.pattern(imputation.data)
final.dat.imp <- mice(
  imputation.data, m=8, maxit=500, method=imputation.method,
  predictorMatrix=imputation.predictors, seed=101
)
summary(final.dat.imp)

final.dat.imp$imp$gdp
final.dat.imp$imp$cpi
final.dat.imp$imp$ghsi
final.dat.imp$imp$igs
final.dat.imp$imp$agrL
final.dat.imp$imp$VAag
final.dat.imp$imp$govexpedu
final.dat.imp$imp$stja

final.dat.compl.list <- lapply(seq_len(8), function(i) {
  completed <- final.dat.raw
  completed[imputation.keep, ] <- complete(final.dat.imp, i)[names(final.dat.raw)]
  completed
})
# Downstream code expects one completed data set.  Preserve all draws above for
# pooled sensitivity analyses rather than replacing missing values by a mean.
final.dat.compl.mn <- final.dat.compl.list[[1]]

gdpsc <- log10(scale((final.dat.compl.mn$gdp), scale=T, center=F))
cpisc <- scale(final.dat.compl.mn$cpi, scale=T, center=F)
ghsisc <- scale((final.dat.compl.mn$ghsi), scale=T, center=F)
igssc <- log10(scale(final.dat.compl.mn$igs, scale=T, center=F))
agrLsc <- scale(logit(final.dat.compl.mn$agrL/100), scale=T, center=F)
VAagsc <- scale(logit(final.dat.compl.mn$VAag/100), scale=T, center=F)
govexpedusc <- scale(logit(final.dat.compl.mn$govexpedu/100), scale=T, center=F)
stjasc <- log10(scale((final.dat.compl.mn$stja), scale=T, center=F))
ratiosc <- log10(scale(ratio.full$dam.mgm, center=F, scale=T))
rsc <- scale(ratio.full$ratio.r, center=F, scale=T)
pdamsc <- scale(logit(ratio.full$pdam), center=F, scale=T)

head(ratio.full)
dim(ratio.full)

final.dat.imputed <- na.omit(data.frame(ratio.full$cntry.code, ratio.full$cont, ratio.full$region, ratiosc,
                                gdpsc, cpisc, ghsisc, igssc, agrLsc, VAagsc, govexpedusc, stjasc))
colnames(final.dat.imputed) <- c("cntry.code", "cont", "reg", "ratio", "gdp", "cpi", "ghsi", "igs", "agrL", "VAag", "govexpedu", "stja")
head(final.dat.imputed)
dim(final.dat.imputed)

final.dat.imputed.orig.scale <- na.omit(data.frame(cntry.code=ratio.full$cntry.code, cont=ratio.full$cont,
                                                   reg=ratio.full$region, ratio=ratio.full$dam.mgm,
                                                   gdp=final.dat.compl.mn$gdp, cpi=final.dat.compl.mn$cpi,
                                                   ghsi=final.dat.compl.mn$ghsi, igs=final.dat.compl.mn$igs,
                                                   agrL=final.dat.compl.mn$agrL, VAag=final.dat.compl.mn$VAag,
                                                   govexpedu=final.dat.compl.mn$govexpedu, stja=final.dat.compl.mn$stja))
head(final.dat.imputed.orig.scale)
colnames(final.dat.imputed.orig.scale)
hist(final.dat.imputed.orig.scale$ratio)
hist(log10(final.dat.imputed.orig.scale$ratio))

hist(sqrt(final.dat.imputed.orig.scale$gdp))
hist(sqrt(final.dat.imputed.orig.scale$cpi))
hist(sqrt(final.dat.imputed.orig.scale$ghsi))
hist(log10(final.dat.imputed.orig.scale$igs))
hist(logit(final.dat.imputed.orig.scale$agrL/100))
hist(logit(final.dat.imputed.orig.scale$VAag/100))
hist((final.dat.imputed.orig.scale$govexpedu))
hist(log10(final.dat.imputed.orig.scale$stja))

final.dat.imputed.orig.scale$ratio.l <- log10(final.dat.imputed.orig.scale$ratio)
final.dat.imputed.orig.scale$gdp.sqr <- sqrt(final.dat.imputed.orig.scale$gdp)
final.dat.imputed.orig.scale$cpi.sqr <- sqrt(final.dat.imputed.orig.scale$cpi)
final.dat.imputed.orig.scale$ghsi.sqr <- sqrt(final.dat.imputed.orig.scale$ghsi)
final.dat.imputed.orig.scale$igs.l <- log10(final.dat.imputed.orig.scale$igs)
final.dat.imputed.orig.scale$agrL.l <- logit(final.dat.imputed.orig.scale$agrL/100)
final.dat.imputed.orig.scale$VAag.l <- logit(final.dat.imputed.orig.scale$VAag/100)
final.dat.imputed.orig.scale$stja.l <- log10(final.dat.imputed.orig.scale$stja)

final.datPdam.imputed <- na.omit(data.frame(ratio.full$cntry.code, ratio.full$cont, ratio.full$region, pdamsc,
                                        gdpsc, cpisc, ghsisc, igssc, agrLsc, VAagsc, govexpedusc, stjasc))
colnames(final.datPdam.imputed) <- c("cntry.code", "cont", "reg", "pdam", "gdp", "cpi", "ghsi", "igs", "agrL", "VAag", "govexpedu", "stja")
head(final.datPdam.imputed)
dim(final.datPdam.imputed)

final.datPdam.imputed.orig.scale <- na.omit(data.frame(cntry.code=ratio.full$cntry.code, cont=ratio.full$cont,
                                                   reg=ratio.full$region, pdam=ratio.full$pdam,
                                                   gdp=final.dat.compl.mn$gdp, cpi=final.dat.compl.mn$cpi,
                                                   ghsi=final.dat.compl.mn$ghsi, igs=final.dat.compl.mn$igs,
                                                   agrL=final.dat.compl.mn$agrL, VAag=final.dat.compl.mn$VAag,
                                                   govexpedu=final.dat.compl.mn$govexpedu, stja=final.dat.compl.mn$stja))
head(final.datPdam.imputed.orig.scale)
colnames(final.datPdam.imputed.orig.scale)
hist(final.datPdam.imputed.orig.scale$pdam)
hist(log10(final.datPdam.imputed.orig.scale$pdam))
hist(logit(final.datPdam.imputed.orig.scale$pdam))

final.datPdam.imputed.orig.scale$pdam.l <- logit(final.datPdam.imputed.orig.scale$pdam)
final.datPdam.imputed.orig.scale$gdp.sqr <- sqrt(final.datPdam.imputed.orig.scale$gdp)
final.datPdam.imputed.orig.scale$cpi.sqr <- sqrt(final.datPdam.imputed.orig.scale$cpi)
final.datPdam.imputed.orig.scale$ghsi.sqr <- sqrt(final.datPdam.imputed.orig.scale$ghsi)
final.datPdam.imputed.orig.scale$igs.l <- log10(final.datPdam.imputed.orig.scale$igs)
final.datPdam.imputed.orig.scale$agrL.l <- logit(final.datPdam.imputed.orig.scale$agrL/100)
final.datPdam.imputed.orig.scale$VAag.l <- logit(final.datPdam.imputed.orig.scale$VAag/100)
final.datPdam.imputed.orig.scale$stja.l <- log10(final.datPdam.imputed.orig.scale$stja)
colnames(final.datPdam.imputed.orig.scale)
head(final.datPdam.imputed.orig.scale)


final.datr.imputed <- na.omit(data.frame(ratio.full$cntry.code, ratio.full$cont, ratio.full$region, rsc,
                                        gdpsc, cpisc, ghsisc, igssc, agrLsc, VAagsc, govexpedusc, stjasc))
colnames(final.datr.imputed) <- c("cntry.code", "cont", "reg", "r", "gdp", "cpi", "ghsi", "igs", "agrL", "VAag", "govexpedu", "stja")
head(final.datr.imputed)
dim(final.datr.imputed)

final.datr.imputed.orig.scale <- na.omit(data.frame(cntry.code=ratio.full$cntry.code, cont=ratio.full$cont,
                                                       reg=ratio.full$region, r=ratio.full$ratio.r,
                                                       gdp=final.dat.compl.mn$gdp, cpi=final.dat.compl.mn$cpi,
                                                       ghsi=final.dat.compl.mn$ghsi, igs=final.dat.compl.mn$igs,
                                                       agrL=final.dat.compl.mn$agrL, VAag=final.dat.compl.mn$VAag,
                                                       govexpedu=final.dat.compl.mn$govexpedu, stja=final.dat.compl.mn$stja))
head(final.datr.imputed.orig.scale)
colnames(final.datr.imputed.orig.scale)
hist(final.datr.imputed.orig.scale$r)

final.datr.imputed.orig.scale$gdp.sqr <- sqrt(final.datr.imputed.orig.scale$gdp)
final.datr.imputed.orig.scale$cpi.sqr <- sqrt(final.datr.imputed.orig.scale$cpi)
final.datr.imputed.orig.scale$ghsi.sqr <- sqrt(final.datr.imputed.orig.scale$ghsi)
final.datr.imputed.orig.scale$igs.l <- log10(final.datr.imputed.orig.scale$igs)
final.datr.imputed.orig.scale$agrL.l <- logit(final.datr.imputed.orig.scale$agrL/100)
final.datr.imputed.orig.scale$VAag.l <- logit(final.datr.imputed.orig.scale$VAag/100)
final.datr.imputed.orig.scale$stja.l <- log10(final.datr.imputed.orig.scale$stja)
colnames(final.datr.imputed.orig.scale)
head(final.datr.imputed.orig.scale)




## correlation matrix
cor.dat.ratio <- data.frame(final.dat.imputed$gdp, final.dat.imputed$cpi, final.dat.imputed$ghsi, final.dat.imputed$igs, final.dat.imputed$agrL,
                            final.dat.imputed$VAag, final.dat.imputed$govexpedu, final.dat.imputed$stja)
colnames(cor.dat.ratio) <- c("GDP","CPI","GHSI","IGS","AGRL","VAPP","EDU","STJA")
cormat.ratio <- cor(na.omit(cor.dat.ratio), method="kendall")
cormat.ratio[lower.tri(cormat.ratio)] <- NA
cormat.ratio
max(abs(cormat.ratio[cormat.ratio < 1]), na.rm=T)
median(abs(cormat.ratio[cormat.ratio < 1]), na.rm=T)
hist(abs(cormat.ratio[cormat.ratio < 1]), main="")

colNam.vec <- colnames(cormat.ratio) # column names
rowNam.vec <- rownames(cormat.ratio) # row names
samp.exp <- expand.grid(X=colNam.vec, Y=rowNam.vec) # expand to all possible row x column pairs
samp.exp$cor <- 0 # set correlation column

# loop through expanded grid to fill with correlation matrix values
for (i in 1:dim(samp.exp)[1]) {
  x.var <- as.character(samp.exp[i,1])
  y.var <- as.character(samp.exp[i,2])
  samp.exp$cor[i] <- cormat.ratio[which(colnames(cormat.ratio) == x.var), which(rownames(cormat.ratio) == y.var)]
}

# plot heatmap
ggplot(samp.exp, aes(X, Y, fill = cor)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", limits=c(-1,1)) +
  labs(x = "", y = "")



#############################
## response: dam:mgm ratio ##
#############################

## deterministic BRT
## variable selection
colnames(final.dat.imputed.orig.scale)
pred.gdp.sub <- which(colnames(final.dat.imputed.orig.scale)=='gdp.sqr')
pred.cpi.sub <- which(colnames(final.dat.imputed.orig.scale)=='cpi.sqr')
pred.ghsi.sub <- which(colnames(final.dat.imputed.orig.scale)=='ghsi.sqr')
pred.igs.sub <- which(colnames(final.dat.imputed.orig.scale)=='igs.l')
pred.agrL.sub <- which(colnames(final.dat.imputed.orig.scale)=='agrL.l')
pred.VAag.sub <- which(colnames(final.dat.imputed.orig.scale)=='VAag.l')
pred.govexpedu.sub <- which(colnames(final.dat.imputed.orig.scale)=='govexpedu')
pred.stja.sub <- which(colnames(final.dat.imputed.orig.scale)=='stja.l')
predictors.sub <- c(pred.gdp.sub,pred.cpi.sub,pred.ghsi.sub,pred.igs.sub,pred.agrL.sub,
                    pred.VAag.sub,pred.govexpedu.sub,pred.stja.sub)
predictors.sub
colnames(final.dat.imputed.orig.scale[,predictors.sub])
resp.sub <- which(colnames(final.dat.imputed.orig.scale)=='ratio.l')
resp.sub
colnames(final.dat.imputed.orig.scale[,c(1,resp.sub)])[2]

brt.fit <- adaptive.gbm.step(final.dat.imputed.orig.scale, gbm.x = attr(final.dat.imputed.orig.scale, "names")[predictors.sub],
                              gbm.y = attr(final.dat.imputed.orig.scale, "names")[resp.sub], 
                             family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.00001, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit)
D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit)
gbm.plot.fits(brt.fit)

brt.CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
brt.CV.cor
brt.CV.cor.se <- 100 * brt.fit$cv.statistics$correlation.se
brt.CV.cor.se
print(c(brt.CV.cor, brt.CV.cor.se))



######################
# resampled BRT loop #
######################
biter <- 1000
eq.sp.points <- 100

## variable selection
colnames(final.dat.imputed.orig.scale)
col_names <- colnames(final.dat.imputed.orig.scale)
pred.gdp.sub <- which(colnames(final.dat.imputed.orig.scale)=='gdp.sqr')
pred.cpi.sub <- which(colnames(final.dat.imputed.orig.scale)=='cpi.sqr')
pred.ghsi.sub <- which(colnames(final.dat.imputed.orig.scale)=='ghsi.sqr')
pred.igs.sub <- which(colnames(final.dat.imputed.orig.scale)=='igs.l')
pred.agrL.sub <- which(colnames(final.dat.imputed.orig.scale)=='agrL.l')
pred.VAag.sub <- which(colnames(final.dat.imputed.orig.scale)=='VAag.l')
pred.govexpedu.sub <- which(colnames(final.dat.imputed.orig.scale)=='govexpedu')
pred.stja.sub <- which(colnames(final.dat.imputed.orig.scale)=='stja.l')
predictors.sub <- c(pred.gdp.sub,pred.cpi.sub,pred.ghsi.sub,pred.igs.sub,pred.agrL.sub,
                    pred.VAag.sub,pred.govexpedu.sub,pred.stja.sub)
predictors.sub
resp.sub <- which(colnames(final.dat.imputed.orig.scale)=='ratio.l')
resp.sub

# create storage arrays
val.arr <- pred.arr <- array(data = 0, dim = c(eq.sp.points, length(predictors.sub), biter),
                             dimnames=list(paste("x",1:eq.sp.points,sep=""),
                            attr(final.dat.imputed.orig.scale, "names")[predictors.sub], paste("b",1:biter,sep="")))

# create storage vectors
D2.vec <- CV.cor.vec <- CV.cor.se.vec <- 
  GDP.ri <- CPI.ri <- GHSI.ri <- IGS.ri <- AGRL.ri <- VAAG.ri <- GOVEXPEDU.ri <- STJA.ri <- rep(0,biter)

for (b in 1:biter) {
  # resample data among countries
  resamp.sub <- sort(sample(x = 1:dim(final.dat.imputed.orig.scale)[1], size = dim(final.dat.imputed.orig.scale)[1], replace=TRUE))
  dat.resamp <- final.dat.imputed.orig.scale[resamp.sub,]
  
  # boosted regression tree
  brt.fit <- adaptive.gbm.step(dat.resamp, gbm.x = attr(dat.resamp, "names")[predictors.sub],
                               gbm.y = attr(dat.resamp, "names")[resp.sub], family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.001, bag.fraction=0.75, tree.complexity = 2, silent=T, tolerance.method = "auto")
  summ.fit <- summary(brt.fit)

  # variable relative importance
  GDP.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[1]])])
  CPI.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[2]])])
  GHSI.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[3]])])
  IGS.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[4]])])
  AGRL.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[5]])])
  VAAG.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[6]])])
  GOVEXPEDU.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[7]])])
  STJA.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[8]])])
  
  D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / 
    brt.fit$cv.statistics$deviance.mean
  D2.vec[b] <- D2
  CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
  CV.cor.vec[b] <- CV.cor
  CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
  CV.cor.se.vec[b] <- CV.cor.se

  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=length(predictors.sub))
  ## output average predictions
  for (p in 1:length(predictors.sub)) {
    RESP.val[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,2]
  }
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.fit$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.fit$var.names
  
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  print(b)
  
} # end b

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- 5
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] > (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  }
  print(k)
}

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)

val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median)

## plot
# order: GDP,CPI,GHSI,IGS,AGRL,VAAG,GOVEXPEDU,STJA
par(mfrow=c(2,4)) 

plot(val.med[,1],pred.med[,1],type="l",ylim=c(min(pred.lo[,1]),max(pred.up[,1])), lwd=2,
     ylab="(←lower) dam:mgm (higher→)", xlab="(←poorer) GDPpc (richer→)" )
lines(val.med[,1], pred.lo[,1], type="l", lty=2, col="red")
lines(val.med[,1], pred.up[,1], type="l", lty=2, col="red")

plot(val.med[,2],pred.med[,2],type="l",ylim=c(min(pred.lo[,2]),max(pred.up[,2])), lwd=2,
     ylab="(←lower) dam:mgm (higher→)", xlab="(←higher) corruption (lower→)" )
lines(val.med[,2], pred.lo[,2], type="l", lty=2, col="red")
lines(val.med[,2], pred.up[,2], type="l", lty=2, col="red")

plot(val.med[,3],pred.med[,3],type="l",ylim=c(min(pred.lo[,3]),max(pred.up[,3])), lwd=2,
     ylab="(←lower) dam:mgm (higher→)", xlab="(←lower) global health security (higher→)" )
lines(val.med[,3], pred.lo[,3], type="l", lty=3, col="red")
lines(val.med[,3], pred.up[,3], type="l", lty=3, col="red")

plot(val.med[,4],pred.med[,4],type="l",ylim=c(min(pred.lo[,4]),max(pred.up[,4])), lwd=2,
     ylab="(←lower) dam:mgm (higher→)", xlab="(←less) goods & services imports (more→)" )
lines(val.med[,4], pred.lo[,4], type="l", lty=3, col="red")
lines(val.med[,4], pred.up[,4], type="l", lty=3, col="red")

plot(val.med[,5],pred.med[,5],type="l",ylim=c(min(pred.lo[,5]),max(pred.up[,5])), lwd=2,
     ylab="(←lower) dam:mgm (higher→)", xlab="(←less) agricultural land (more→)" )
lines(val.med[,5], pred.lo[,5], type="l", lty=2, col="red")
lines(val.med[,5], pred.up[,5], type="l", lty=2, col="red")

plot(val.med[,6],pred.med[,6],type="l",ylim=c(min(pred.lo[,6]),max(pred.up[,6])), lwd=2,
     ylab="(←lower) dam:mgm (higher→)", xlab="(←less) value added primary production (more→)" )
lines(val.med[,6], pred.lo[,6], type="l", lty=2, col="red")
lines(val.med[,6], pred.up[,6], type="l", lty=2, col="red")

plot(val.med[,7],pred.med[,7],type="l",ylim=c(min(pred.lo[,7]),max(pred.up[,7])), lwd=2,
     ylab="(←lower) dam:mgm (higher→)", xlab="(←less) education investment (more→)" )
lines(val.med[,7], pred.lo[,7], type="l", lty=2, col="red")
lines(val.med[,7], pred.up[,7], type="l", lty=2, col="red")

plot(val.med[,8],pred.med[,8],type="l",ylim=c(min(pred.lo[,8]),max(pred.up[,8])), lwd=2,
     ylab="(←lower) dam:mgm (higher→)", xlab="(←fewer) scientific & tech journal articles (more→)" )
lines(val.med[,8], pred.lo[,8], type="l", lty=2, col="red")
lines(val.med[,8], pred.up[,8], type="l", lty=2, col="red")

par(mfrow=c(1,1)) 

# kappa method for output vectors
D2.update <- D2.vec[1:biter]
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# order: GDP,CPI,GHSI,IGS,AGRL,VAAG,GOVEXPEDU,STJA
GDP.ri.update <- GDP.ri[1:biter]
CPI.ri.update <- CPI.ri[1:biter]
GHSI.ri.update <- GHSI.ri[1:biter]
IGS.ri.update <- IGS.ri[1:biter]
AGRL.ri.update <- AGRL.ri[1:biter]
VAAG.ri.update <- VAAG.ri[1:biter]
GOVEXPEDU.ri.update <- GOVEXPEDU.ri[1:biter]
STJA.ri.update <- STJA.ri[1:biter]

for (k in 1:kappa.n) {
  D2.mean <- mean(D2.update, na.rm=T); D2.sd <- sd(D2.update, na.rm=T)
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  GDP.mean <- mean(GDP.ri.update, na.rm=T); GDP.sd <- sd(GDP.ri.update, na.rm=T)
  CPI.mean <- mean(CPI.ri.update, na.rm=T); CPI.sd <- sd(CPI.ri.update, na.rm=T)
  GHSI.mean <- mean(GHSI.ri.update, na.rm=T); GHSI.sd <- sd(GHSI.ri.update, na.rm=T)
  IGS.mean <- mean(IGS.ri.update, na.rm=T); IGS.sd <- sd(IGS.ri.update, na.rm=T)
  AGRL.mean <- mean(AGRL.ri.update, na.rm=T); AGRL.sd <- sd(AGRL.ri.update, na.rm=T)
  VAAG.mean <- mean(VAAG.ri.update, na.rm=T); VAAG.sd <- sd(VAAG.ri.update, na.rm=T)
  GOVEXPEDU.mean <- mean(GOVEXPEDU.ri.update, na.rm=T); GOVEXPEDU.sd <- sd(GOVEXPEDU.ri.update, na.rm=T)
  STJA.mean <- mean(STJA.ri.update, na.rm=T); STJA.sd <- sd(STJA.ri.update, na.rm=T)

  for (u in 1:biter) {
    D2.update[u] <- ifelse((D2.update[u] < (D2.mean-kappa*D2.sd) | D2.update[u] > (D2.mean+kappa*D2.sd)), NA, D2.update[u])
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] > (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] > (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    
    # order: GDP,CPI,GHSI,IGS,AGRL,VAAG,GOVEXPEDU,STJA
    GDP.ri.update[u] <- ifelse((GDP.ri.update[u] < (GDP.mean-kappa*GDP.sd) | GDP.ri.update[u] > (GDP.mean+kappa*GDP.sd)), NA, GDP.ri.update[u])
    CPI.ri.update[u] <- ifelse((CPI.ri.update[u] < (CPI.mean-kappa*CPI.sd) | CPI.ri.update[u] > (CPI.mean+kappa*CPI.sd)), NA, CPI.ri.update[u])
    GHSI.ri.update[u] <- ifelse((GHSI.ri.update[u] < (GHSI.mean-kappa*GHSI.sd) | GHSI.ri.update[u] > (GHSI.mean+kappa*GHSI.sd)), NA, GHSI.ri.update[u])
    IGS.ri.update[u] <- ifelse((IGS.ri.update[u] < (IGS.mean-kappa*IGS.sd) | IGS.ri.update[u] > (IGS.mean+kappa*IGS.sd)), NA, IGS.ri.update[u])
    AGRL.ri.update[u] <- ifelse((AGRL.ri.update[u] < (AGRL.mean-kappa*AGRL.sd) | AGRL.ri.update[u] > (AGRL.mean+kappa*AGRL.sd)), NA, AGRL.ri.update[u])
    VAAG.ri.update[u] <- ifelse((VAAG.ri.update[u] < (VAAG.mean-kappa*VAAG.sd) | VAAG.ri.update[u] > (VAAG.mean+kappa*VAAG.sd)), NA, VAAG.ri.update[u])
    GOVEXPEDU.ri.update[u] <- ifelse((GOVEXPEDU.ri.update[u] < (GOVEXPEDU.mean-kappa*GOVEXPEDU.sd) | GOVEXPEDU.ri.update[u] > (GOVEXPEDU.mean+kappa*GOVEXPEDU.sd)), NA, GOVEXPEDU.ri.update[u])
    STJA.ri.update[u] <- ifelse((STJA.ri.update[u] < (STJA.mean-kappa*STJA.sd) | STJA.ri.update[u] > (STJA.mean+kappa*STJA.sd)), NA, STJA.ri.update[u])
  }
  print(k)
}

D2.med <- median(D2.update, na.rm=TRUE)
D2.lo <- quantile(D2.update, probs=0.025, na.rm=TRUE)
D2.up <- quantile(D2.update, probs=0.975, na.rm=TRUE)
print(c(D2.lo,D2.med,D2.up))

CV.cor.med <- median(CV.cor.update, na.rm=TRUE)
CV.cor.lo <- quantile(CV.cor.update, probs=0.025, na.rm=TRUE)
CV.cor.up <- quantile(CV.cor.update, probs=0.975, na.rm=TRUE)
print(c(CV.cor.lo,CV.cor.med,CV.cor.up))

# order: GDP,CPI,GHSI,IGS,AGRL,VAAG,GOVEXPEDU,STJA
GDP.ri.lo <- quantile(GDP.ri.update, probs=0.025, na.rm=TRUE)
GDP.ri.med <- median(GDP.ri.update, na.rm=TRUE)
GDP.ri.up <- quantile(GDP.ri.update, probs=0.975, na.rm=TRUE)

CPI.ri.lo <- quantile(CPI.ri.update, probs=0.025, na.rm=TRUE)
CPI.ri.med <- median(CPI.ri.update, na.rm=TRUE)
CPI.ri.up <- quantile(CPI.ri.update, probs=0.975, na.rm=TRUE)

GHSI.ri.lo <- quantile(GHSI.ri.update, probs=0.025, na.rm=TRUE)
GHSI.ri.med <- median(GHSI.ri.update, na.rm=TRUE)
GHSI.ri.up <- quantile(GHSI.ri.update, probs=0.975, na.rm=TRUE)

IGS.ri.lo <- quantile(IGS.ri.update, probs=0.025, na.rm=TRUE)
IGS.ri.med <- median(IGS.ri.update, na.rm=TRUE)
IGS.ri.up <- quantile(IGS.ri.update, probs=0.975, na.rm=TRUE)

AGRL.ri.lo <- quantile(AGRL.ri.update, probs=0.025, na.rm=TRUE)
AGRL.ri.med <- median(AGRL.ri.update, na.rm=TRUE)
AGRL.ri.up <- quantile(AGRL.ri.update, probs=0.975, na.rm=TRUE)

VAAG.ri.lo <- quantile(VAAG.ri.update, probs=0.025, na.rm=TRUE)
VAAG.ri.med <- median(VAAG.ri.update, na.rm=TRUE)
VAAG.ri.up <- quantile(VAAG.ri.update, probs=0.975, na.rm=TRUE)

GOVEXPEDU.ri.lo <- quantile(GOVEXPEDU.ri.update, probs=0.025, na.rm=TRUE)
GOVEXPEDU.ri.med <- median(GOVEXPEDU.ri.update, na.rm=TRUE)
GOVEXPEDU.ri.up <- quantile(GOVEXPEDU.ri.update, probs=0.975, na.rm=TRUE)

STJA.ri.lo <- quantile(STJA.ri.update, probs=0.025, na.rm=TRUE)
STJA.ri.med <- median(STJA.ri.update, na.rm=TRUE)
STJA.ri.up <- quantile(STJA.ri.update, probs=0.975, na.rm=TRUE)

# order: GDP,CPI,GHSI,IGS,AGRL,VAAG,GOVEXPEDU,STJA
ri.lo <- c(GDP.ri.lo,CPI.ri.lo,GHSI.ri.lo,IGS.ri.lo,AGRL.ri.lo,VAAG.ri.lo,GOVEXPEDU.ri.lo,STJA.ri.lo)
ri.med <- c(GDP.ri.med,CPI.ri.med,GHSI.ri.med,IGS.ri.med,AGRL.ri.med,VAAG.ri.med,GOVEXPEDU.ri.med,STJA.ri.med)
ri.up <- c(GDP.ri.up,CPI.ri.up,GHSI.ri.up,IGS.ri.up,AGRL.ri.up,VAAG.ri.up,GOVEXPEDU.ri.up,STJA.ri.up)

ri.out <- as.data.frame(cbind(ri.lo,ri.med,ri.up))
colnames(ri.out) <- c("ri.lo","ri.med","ri.up")
rownames(ri.out) <- attr(final.dat.imputed.orig.scale, "names")[predictors.sub]
ri.sort <- ri.out[order(ri.out[,2],decreasing=T),1:3]
ri.sort

ratio.partial.pred.med <- pred.med
ratio.partial.pred.up <- pred.up
ratio.partial.pred.lo <- pred.lo
ratio.partial.val.med <- val.med

# write outputs
dir.path <- paste(getwd(),"/out/",sep="")
for (o in 1:length(predictors.sub)) {
  var.name <- colnames(final.dat.imputed.orig.scale)[predictors.sub[o]]
  var.name.abbr <- toupper(sub("\\..*", "", var.name, ignore.case = T))
  var.title <- paste("ratio.partial.",var.name.abbr,sep="")
  assign(var.title, 
         data.frame(val=ratio.partial.val.med[,o],med=as.numeric(ratio.partial.pred.med[,o]), 
                    up=as.numeric(ratio.partial.pred.up[,o]), lo=as.numeric(ratio.partial.pred.lo[,o])))
  get(var.title)
  write.table(get(var.title),paste(dir.path,paste(var.title,".csv",sep=""),sep=""), sep=",", row.names=F, col.names=T)
}

var.names <- rownames(ri.sort)
var.names.abbr <- toupper(sub("\\..*", "", var.names, ignore.case = T))
rel.infl.out <- data.frame(var=var.names.abbr, med=ri.sort$ri.med, up=ri.sort$ri.up, lo=ri.sort$ri.lo)
write.table(rel.infl.out,paste(dir.path,"ratio.rel.infl.csv", sep=""), sep=",", row.names = F, col.names = T)



##########################
## phase-based approach ##
##########################
## response: dam:mgm ratio
## PHASE 1
##########
colnames(final.dat.imputed.orig.scale)
pred.gdp.sub <- which(colnames(final.dat.imputed.orig.scale)=='gdp.sqr')
pred.cpi.sub <- which(colnames(final.dat.imputed.orig.scale)=='cpi.sqr')
pred.ghsi.sub <- which(colnames(final.dat.imputed.orig.scale)=='ghsi.sqr')
pred.igs.sub <- which(colnames(final.dat.imputed.orig.scale)=='igs.l')
predictors.sub <- c(pred.gdp.sub,pred.cpi.sub,pred.ghsi.sub,pred.igs.sub)
resp.sub <- which(colnames(final.dat.imputed.orig.scale)=='ratio.l')

brt.fit1 <- adaptive.gbm.step(final.dat.imputed.orig.scale, gbm.x = attr(final.dat.imputed.orig.scale, "names")[predictors.sub],
                              gbm.y = attr(final.dat.imputed.orig.scale, "names")[resp.sub], family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.00001, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit1)
D2 <- 100 * (brt.fit1$cv.statistics$deviance.mean - brt.fit1$self.statistics$mean.resid) / brt.fit1$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit1)
gbm.plot.fits(brt.fit1)

brt1.CV.cor <- 100 * brt.fit1$cv.statistics$correlation.mean
brt1.CV.cor
brt1.CV.cor.se <- 100 * brt.fit1$cv.statistics$correlation.se
brt1.CV.cor.se
print(c(brt1.CV.cor, brt1.CV.cor.se))


# GLMM
# recode continent to greater region to increase per-category sample sizes
table(final.dat.imputed.orig.scale$cont)
final.dat.imputed.orig.scale$reg2 <- ifelse(final.dat.imputed.orig.scale$cont == "NAM" | final.dat.imputed.orig.scale$cont == "CAR", "NAMCAR", final.dat.imputed.orig.scale$cont)
final.dat.imputed.orig.scale$reg2 <- ifelse(final.dat.imputed.orig.scale$cont == "EUR" | final.dat.imputed.orig.scale$cont == "ME", "EURME", final.dat.imputed.orig.scale$reg2)
final.dat.imputed.orig.scale$reg2 <- ifelse(final.dat.imputed.orig.scale$cont == "ASIA" | final.dat.imputed.orig.scale$cont == "OC", "ASIAOC", final.dat.imputed.orig.scale$reg2)
table(final.dat.imputed.orig.scale$reg2)

# model set
colnames(final.dat.imputed.orig.scale)
vars <- c("gdp.sqr","cpi.sqr","ghsi.sqr","igs.l")    
vars4comb <- paste("ratio.l~", apply(combn(vars,4),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars3comb <- paste("ratio.l~", apply(combn(vars,3),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars2comb <- paste("ratio.l~", apply(combn(vars,2),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars1comb <- paste("ratio.l~", vars, "+(1|reg2)", sep="")
interconly <- "ratio.l~1+(1|reg2)"
mod.vec <- c(vars4comb,vars3comb,vars2comb,vars1comb,interconly)
mod.vec  

length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.dat.imputed.orig.scale, na.action=na.omit, REML=FALSE, verbose=F)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- AICc(fit)
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r2_nakagawa(fit)$R2_marginal # marginal R-squared
  Rc[i] <- 100*r2_nakagawa(fit)$R2_conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable1 <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable1) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable1) <- as.character(mod.vec)
summary.table1 <- sumtable1[order(sumtable1[,5],decreasing=F),]
summary.table1

fitsat <- lmer(as.formula(mod.vec[1]), data=final.dat.imputed.orig.scale, na.action=na.omit, REML=FALSE)
check_model(fitsat, detrend=T)
plot_model(fitsat, show.values=T, vline.color = "purple")
vif(fitsat)

## gdp collinear; remove
# model set
colnames(final.dat.imputed.orig.scale)
vars <- c("cpi.sqr","ghsi.sqr","igs.l")    
vars3comb <- paste("ratio.l~", apply(combn(vars,3),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars2comb <- paste("ratio.l~", apply(combn(vars,2),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars1comb <- paste("ratio.l~", vars, "+(1|reg2)", sep="")
interconly <- "ratio.l~1+(1|reg2)"
mod.vec <- c(vars3comb,vars2comb,vars1comb,interconly)
mod.vec  

length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.dat.imputed.orig.scale, na.action=na.omit, REML=FALSE, verbose=F)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- AICc(fit)
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r2_nakagawa(fit)$R2_marginal # marginal R-squared
  Rc[i] <- 100*r2_nakagawa(fit)$R2_conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable1 <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable1) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable1) <- as.character(mod.vec)
summary.table1 <- sumtable1[order(sumtable1[,5],decreasing=F),]
summary.table1

fitsat <- lmer(as.formula(mod.vec[1]), data=final.dat.imputed.orig.scale, na.action=na.omit, REML=FALSE)
check_model(fitsat, detrend=T)
plot_model(fitsat, show.values=T, vline.color = "purple")
vif(fitsat)


## PHASE 2 - primary production
# agrL, VAag
colnames(final.dat.imputed.orig.scale)
pred.agrl.sub <- which(colnames(final.dat.imputed.orig.scale)=='agrL.l')
pred.vaag.sub <- which(colnames(final.dat.imputed.orig.scale)=='VAag.l')
predictors.sub <- c(pred.agrl.sub,pred.vaag.sub)
resp.sub <- which(colnames(final.dat.imputed.orig.scale)=='ratio.l')

brt.fit1 <- adaptive.gbm.step(final.dat.imputed.orig.scale, gbm.x = attr(final.dat.imputed.orig.scale, "names")[predictors.sub],
                              gbm.y = attr(final.dat.imputed.orig.scale, "names")[resp.sub], family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.00001, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit1)
D2 <- 100 * (brt.fit1$cv.statistics$deviance.mean - brt.fit1$self.statistics$mean.resid) / brt.fit1$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit1)
gbm.plot.fits(brt.fit1)

brt1.CV.cor <- 100 * brt.fit1$cv.statistics$correlation.mean
brt1.CV.cor
brt1.CV.cor.se <- 100 * brt.fit1$cv.statistics$correlation.se
brt1.CV.cor.se
print(c(brt1.CV.cor, brt1.CV.cor.se))

# model set
vars <- c("agrL.l","VAag.l")    
vars2comb <- paste("ratio.l~", apply(combn(vars,2),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars1comb <- paste("ratio.l~", vars, "+(1|reg2)", sep="")
interconly <- "ratio.l~1+(1|reg2)"
mod.vec <- c(vars2comb,vars1comb,interconly)
mod.vec  

length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.dat.imputed.orig.scale, na.action=na.omit, REML=FALSE)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- AICc(fit)
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r2_nakagawa(fit)$R2_marginal # marginal R-squared
  Rc[i] <- 100*r2_nakagawa(fit)$R2_conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable2 <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable2) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable2) <- as.character(mod.vec)
summary.table2 <- sumtable2[order(sumtable2[,5],decreasing=F),]
summary.table2

fitsat <- lmer(as.formula(mod.vec[1]), data=final.dat.imputed.orig.scale, na.action=na.omit, REML=FALSE)
check_model(fitsat, detrend=T)
plot_model(fitsat, show.values=T, vline.color = "purple")


## PHASE 3 - education/research
# govexpedu, stja
colnames(final.dat.imputed.orig.scale)
pred.govexpedu.sub <- which(colnames(final.dat.imputed.orig.scale)=='govexpedu')
pred.stja.l.sub <- which(colnames(final.dat.imputed.orig.scale)=='stja.l')
predictors.sub <- c(pred.govexpedu.sub,pred.stja.l.sub)
resp.sub <- which(colnames(final.dat.imputed.orig.scale)=='ratio.l')

brt.fit1 <- adaptive.gbm.step(final.dat.imputed.orig.scale, gbm.x = attr(final.dat.imputed.orig.scale, "names")[predictors.sub],
                              gbm.y = attr(final.dat.imputed.orig.scale, "names")[resp.sub], family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.00001, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit1)
D2 <- 100 * (brt.fit1$cv.statistics$deviance.mean - brt.fit1$self.statistics$mean.resid) / brt.fit1$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit1)
gbm.plot.fits(brt.fit1)

brt1.CV.cor <- 100 * brt.fit1$cv.statistics$correlation.mean
brt1.CV.cor
brt1.CV.cor.se <- 100 * brt.fit1$cv.statistics$correlation.se
brt1.CV.cor.se
print(c(brt1.CV.cor, brt1.CV.cor.se))


# GLMM
# model set
colnames(final.dat.imputed.orig.scale)
vars <- c("govexpedu","stja.l")    
vars2comb <- paste("ratio.l~", apply(combn(vars,2),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars1comb <- paste("ratio.l~", vars, "+(1|reg2)", sep="")
interconly <- "ratio.l~1+(1|reg2)"
mod.vec <- c(vars2comb,vars1comb,interconly)
mod.vec  

length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.dat.imputed.orig.scale, na.action=na.omit, REML=FALSE)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- AIC(fit)
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r2_nakagawa(fit)$R2_marginal # marginal R-squared
  Rc[i] <- 100*r2_nakagawa(fit)$R2_conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable3 <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable3) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable3) <- as.character(mod.vec)
summary.table3 <- sumtable3[order(sumtable3[,5],decreasing=F),]
summary.table3

fitsat <- lmer(as.formula(mod.vec[1]), data=final.dat.imputed.orig.scale, na.action=na.omit, REML=FALSE)
check_model(fitsat, detrend=T)
plot_model(fitsat, show.values=T, vline.color = "purple")

# remove stja
# model set
colnames(final.dat.imputed.orig.scale)
vars <- c("govexpedu")    
vars1comb <- paste("ratio.l~", vars, "+(1|reg2)", sep="")
interconly <- "ratio.l~1+(1|reg2)"
mod.vec <- c(vars1comb,interconly)
mod.vec  

length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.dat.imputed.orig.scale, na.action=na.omit, REML=FALSE)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- AIC(fit)
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r2_nakagawa(fit)$R2_marginal # marginal R-squared
  Rc[i] <- 100*r2_nakagawa(fit)$R2_conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable3 <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable3) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable3) <- as.character(mod.vec)
summary.table3 <- sumtable3[order(sumtable3[,5],decreasing=F),]
summary.table3

fitsat <- lmer(as.formula(mod.vec[1]), data=final.dat.imputed.orig.scale, na.action=na.omit, REML=FALSE)
check_model(fitsat, detrend=T)
plot_model(fitsat, show.values=T, vline.color = "purple")




## COMBINE PHASES
# ghsi, igs, cpi, agrL, govexpedu
# model set
colnames(final.dat.imputed.orig.scale)
vars <- c("ghsi.sqr","igs.l","cpi.sqr","agrL.l","govexpedu")    
vars5comb <- paste("ratio.l~", apply(combn(vars,5),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars4comb <- paste("ratio.l~", apply(combn(vars,4),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars3comb <- paste("ratio.l~", apply(combn(vars,3),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars2comb <- paste("ratio.l~", apply(combn(vars,2),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars1comb <- paste("ratio.l~", vars, "+(1|reg2)", sep="")
interconly <- "ratio.l~1+(1|reg2)"
mod.vec <- c(vars5comb,vars4comb,vars3comb,vars2comb,vars1comb,interconly)
mod.vec  

length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.dat.imputed.orig.scale, na.action=na.omit, REML=FALSE)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- AICc(fit)
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r2_nakagawa(fit)$R2_marginal # marginal R-squared
  Rc[i] <- 100*r2_nakagawa(fit)$R2_conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtableF <- data.frame(mod.num,k.vec,round(LL.vec,3),round(AICc.vec,3),round(dAICc,3),round(wAICc,3),
                        round(BIC.vec,3),round(dBIC,3),round(wBIC,4),round(Rm,1),round(Rc,1))
colnames(sumtableF) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtableF) <- as.character(mod.vec)
summary.tableF <- sumtableF[order(sumtableF[,5],decreasing=F),]
summary.tableF

fitsat <- lmer(as.formula(mod.vec[1]), data=final.dat.imputed.orig.scale, na.action=na.omit, REML=FALSE)
vif(fitsat)
check_model(fitsat, detrend=F)
plot_model(fitsat)
plot_model(fitsat, type="re")


## BRT
# ghsi, igs, cpi, agrL, govexpedu
colnames(final.dat.imputed.orig.scale)
pred.ghsi.sub <- which(colnames((final.dat.imputed.orig.scale))=='ghsi.sqr')
pred.igs.sub <- which(colnames(final.dat.imputed.orig.scale)=='igs.l')
pred.cpi.sub <- which(colnames(final.dat.imputed.orig.scale)=='cpi.sqr')
pred.agrL.sub <- which(colnames(final.dat.imputed.orig.scale)=='agrL.l')
pred.govexpedu.sub <- which(colnames(final.dat.imputed.orig.scale)=='govexpedu')
predictors.sub <- c(pred.ghsi.sub,pred.igs.sub,pred.cpi.sub,pred.agrL.sub,pred.govexpedu.sub)
resp.sub <- which(colnames(final.dat.imputed.orig.scale)=='ratio.l')

brt.fitF <- adaptive.gbm.step(final.dat.imputed.orig.scale, gbm.x = attr(final.dat.imputed.orig.scale, "names")[predictors.sub],
                              gbm.y = attr(final.dat.imputed.orig.scale, "names")[resp.sub], family="gaussian", max.trees=100000, tolerance = 0.00001, learning.rate = 0.00003, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fitF)
D2 <- 100 * (brt.fitF$cv.statistics$deviance.mean - brt.fitF$self.statistics$mean.resid) / brt.fitF$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fitF)
gbm.plot.fits(brt.fitF)

brtF.CV.cor <- 100 * brt.fitF$cv.statistics$correlation.mean
brtF.CV.cor
brtF.CV.cor.se <- 100 * brt.fitF$cv.statistics$correlation.se
brtF.CV.cor.se
print(c(brtF.CV.cor, brtF.CV.cor.se))



# resampled BRT loop
# ghsi, igs, cpi, agrL, govexpedu
biter <- 1000
eq.sp.points <- 100

colnames(final.dat.imputed.orig.scale)
pred.ghsi.sub <- which(colnames((final.dat.imputed.orig.scale))=='ghsi.sqr')
pred.igs.sub <- which(colnames(final.dat.imputed.orig.scale)=='igs.l')
pred.cpi.sub <- which(colnames(final.dat.imputed.orig.scale)=='cpi.sqr')
pred.agrL.sub <- which(colnames(final.dat.imputed.orig.scale)=='agrL.l')
pred.govexpedu.sub <- which(colnames(final.dat.imputed.orig.scale)=='govexpedu')
predictors.sub <- c(pred.ghsi.sub,pred.igs.sub,pred.cpi.sub,pred.agrL.sub,pred.govexpedu.sub)
resp.sub <- which(colnames(final.dat.imputed.orig.scale)=='ratio.l')

# create storage arrays
val.arr <- pred.arr <- array(data = NA, dim = c(eq.sp.points, length(predictors.sub), biter),
                             dimnames=list(paste("x",1:eq.sp.points,sep=""), 
                             attr(final.dat.imputed.orig.scale, "names")[predictors.sub], paste("b",1:biter,sep="")))

# create storage vectors
D2.vec <- CV.cor.vec <- CV.cor.se.vec <- GHSI.ri <- IGS.ri <- CPI.ri <- AGRIL.ri <- GOVEXPEDU.ri <- rep(NA,biter)

for (b in 1:biter) {
  # resample data among countries
  resamp.sub <- sort(sample(x = 1:dim(final.dat.imputed.orig.scale)[1], size = dim(final.dat.imputed.orig.scale)[1], replace=TRUE))
  dat.resamp <- final.dat.imputed.orig.scale[resamp.sub,]
  
  # boosted regression tree
  brt.fit <- adaptive.gbm.step(dat.resamp, gbm.x = attr(dat.resamp, "names")[predictors.sub],
                               gbm.y = attr(dat.resamp, "names")[resp.sub], family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.0001, bag.fraction=0.75, tree.complexity = 2, silent=T, tolerance.method = "auto")
  summ.fit <- summary(brt.fit)
  length(summ.fit[[1]])
  
    # variable relative importance
    GHSI.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[1]])])
    IGS.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[2]])])
    CPI.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[3]])])
    AGRL.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[4]])])
    GOVEXPEDU.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[5]])])
    
    D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
    D2.vec[b] <- D2
    CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
    CV.cor.vec[b] <- CV.cor
    CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
    CV.cor.se.vec[b] <- CV.cor.se
    
    RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=length(predictors.sub))
    ## output average predictions
    for (p in 1:length(predictors.sub)) {
      RESP.val[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,1]
      RESP.pred[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,2]
    }
    RESP.val.dat <- as.data.frame(RESP.val)
    colnames(RESP.val.dat) <- brt.fit$var.names
    RESP.pred.dat <- as.data.frame(RESP.pred)
    colnames(RESP.pred.dat) <- brt.fit$var.names
    
    val.arr[, , b] <- as.matrix(RESP.val.dat)
    pred.arr[, , b] <- as.matrix(RESP.pred.dat)
    
    print(b)
  
} # end b

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- 5
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  }
  print(k)
}

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)

val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

par(mfrow=c(2,3)) 
# igs, agrL, govexpedu
plot(val.med[,1],pred.med[,1],type="l",ylim=c(min(pred.lo[,1]),max(pred.up[,1])),
     lwd=2, ylab="(←lower) dam:mgm (higher→)",  xlab="(←lower) global health security index (higher→)" )
lines(val.med[,1], pred.lo[,1], type="l", lty=2, col="red")
lines(val.med[,1], pred.up[,1], type="l", lty=2, col="red")

plot(val.med[,2],pred.med[,2],type="l",ylim=c(min(pred.lo[,2]),max(pred.up[,2])),
     lwd=2, ylab="(←lower) dam:mgm (higher→)",  xlab="(←lower) per capita imports goods & services (higher→)" )
lines(val.med[,2], pred.lo[,2], type="l", lty=2, col="red")
lines(val.med[,2], pred.up[,2], type="l", lty=2, col="red")

plot(val.med[,3],pred.med[,3],type="l",ylim=c(min(pred.lo[,3]),max(pred.up[,3])),
     lwd=2, ylab="(←lower) dam:mgm (higher→)",  xlab="(←higher) corruption perception index (lower→)" )
lines(val.med[,3], pred.lo[,3], type="l", lty=2, col="red")
lines(val.med[,3], pred.up[,3], type="l", lty=2, col="red")

plot(val.med[,4],pred.med[,4],type="l",ylim=c(min(pred.lo[,4]),max(pred.up[,4])),
     lwd=3, ylab="(←lower) dam:mgm (higher→)", xlab="(←less) agricultural land (more→)")
lines(val.med[,4], pred.lo[,4], type="l", lty=2, col="red")
lines(val.med[,4], pred.up[,4], type="l", lty=2, col="red")

plot(val.med[,5],pred.med[,5],type="l",ylim=c(min(pred.lo[,5]),max(pred.up[,5])),
     lwd=3, ylab="(←lower) dam:mgm (higher→)", xlab="(←less) gov education expenditure (more→)")
lines(val.med[,5], pred.lo[,5], type="l", lty=2, col="red")
lines(val.med[,5], pred.up[,5], type="l", lty=2, col="red")

par(mfrow=c(1,1)) 

# kappa method for output vectors
D2.update <- D2.vec[1:biter]
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# igs, agrL, govexpedu
GHSI.ri.update <- GHSI.ri[1:biter]
IGS.ri.update <- IGS.ri[1:biter]
CPI.ri.update <- CPI.ri[1:biter]
AGRL.ri.update <- AGRL.ri[1:biter]
GOVEXPEDU.ri.update <- GOVEXPEDU.ri[1:biter]

for (k in 1:kappa.n) {
  D2.mean <- mean(D2.update, na.rm=T); D2.sd <- sd(D2.update, na.rm=T)
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  GHSI.mean <- mean(GHSI.ri.update, na.rm=T); GHSI.sd <- sd(GHSI.ri.update, na.rm=T)
  IGS.mean <- mean(IGS.ri.update, na.rm=T); IGS.sd <- sd(IGS.ri.update, na.rm=T)
  CPI.mean <- mean(CPI.ri.update, na.rm=T); CPI.sd <- sd(CPI.ri.update, na.rm=T)
  AGRL.mean <- mean(AGRL.ri.update, na.rm=T); AGRL.sd <- sd(AGRL.ri.update, na.rm=T)
  GOVEXPEDU.mean <- mean(GOVEXPEDU.ri.update, na.rm=T); GOVEXPEDU.sd <- sd(GOVEXPEDU.ri.update, na.rm=T)
  
  for (u in 1:biter) {
    D2.update[u] <- ifelse((D2.update[u] < (D2.mean-kappa*D2.sd) | D2.update[u] > (D2.mean+kappa*D2.sd)), NA, D2.update[u])
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] > (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] > (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    
    GHSI.ri.update[u] <- ifelse((GHSI.ri.update[u] < (GHSI.mean-kappa*GHSI.sd) | GHSI.ri.update[u] > (GHSI.mean+kappa*GHSI.sd)), NA, GHSI.ri.update[u])
    IGS.ri.update[u] <- ifelse((IGS.ri.update[u] < (IGS.mean-kappa*IGS.sd) | IGS.ri.update[u] > (IGS.mean+kappa*IGS.sd)), NA, IGS.ri.update[u])
    CPI.ri.update[u] <- ifelse((CPI.ri.update[u] < (CPI.mean-kappa*CPI.sd) | CPI.ri.update[u] > (CPI.mean+kappa*CPI.sd)), NA, CPI.ri.update[u])
    AGRL.ri.update[u] <- ifelse((AGRL.ri.update[u] < (AGRL.mean-kappa*AGRL.sd) | AGRL.ri.update[u] > (AGRL.mean+kappa*AGRL.sd)), NA, AGRL.ri.update[u])
    GOVEXPEDU.ri.update[u] <- ifelse((GOVEXPEDU.ri.update[u] < (GOVEXPEDU.mean-kappa*GOVEXPEDU.sd) | GOVEXPEDU.ri.update[u] > (GOVEXPEDU.mean+kappa*GOVEXPEDU.sd)), NA, GOVEXPEDU.ri.update[u])
  }
  
  print(k)
}

D2.med <- median(D2.update, na.rm=TRUE)
D2.lo <- quantile(D2.update, probs=0.025, na.rm=TRUE)
D2.up <- quantile(D2.update, probs=0.975, na.rm=TRUE)
print(c(D2.lo,D2.med,D2.up))

CV.cor.med <- median(CV.cor.update, na.rm=TRUE)
CV.cor.lo <- quantile(CV.cor.update, probs=0.025, na.rm=TRUE)
CV.cor.up <- quantile(CV.cor.update, probs=0.975, na.rm=TRUE)
print(c(CV.cor.lo,CV.cor.med,CV.cor.up))

GHSI.ri.lo <- quantile(GHSI.ri.update, probs=0.025, na.rm=TRUE)
GHSI.ri.med <- median(GHSI.ri.update, na.rm=TRUE)
GHSI.ri.up <- quantile(GHSI.ri.update, probs=0.975, na.rm=TRUE)

IGS.ri.lo <- quantile(IGS.ri.update, probs=0.025, na.rm=TRUE)
IGS.ri.med <- median(IGS.ri.update, na.rm=TRUE)
IGS.ri.up <- quantile(IGS.ri.update, probs=0.975, na.rm=TRUE)

CPI.ri.lo <- quantile(CPI.ri.update, probs=0.025, na.rm=TRUE)
CPI.ri.med <- median(CPI.ri.update, na.rm=TRUE)
CPI.ri.up <- quantile(CPI.ri.update, probs=0.975, na.rm=TRUE)

AGRL.ri.lo <- quantile(AGRL.ri.update, probs=0.025, na.rm=TRUE)
AGRL.ri.med <- median(AGRL.ri.update, na.rm=TRUE)
AGRL.ri.up <- quantile(AGRL.ri.update, probs=0.975, na.rm=TRUE)

GOVEXPEDU.ri.lo <- quantile(GOVEXPEDU.ri.update, probs=0.025, na.rm=TRUE)
GOVEXPEDU.ri.med <- median(GOVEXPEDU.ri.update, na.rm=TRUE)
GOVEXPEDU.ri.up <- quantile(GOVEXPEDU.ri.update, probs=0.975, na.rm=TRUE)

ri.lo <- c(GHSI.ri.lo,IGS.ri.lo,CPI.ri.lo,AGRL.ri.lo,GOVEXPEDU.ri.lo)
ri.med <- c(GHSI.ri.med,IGS.ri.med,CPI.ri.med,AGRL.ri.med,GOVEXPEDU.ri.med)
ri.up <- c(GHSI.ri.up,IGS.ri.up,CPI.ri.up,AGRL.ri.up,GOVEXPEDU.ri.up)

ri.out <- as.data.frame(cbind(ri.lo,ri.med,ri.up))
colnames(ri.out) <- c("ri.lo","ri.med","ri.up")
rownames(ri.out) <- attr(final.dat.imputed.orig.scale, "names")[predictors.sub]
ri.sort <- ri.out[order(ri.out[,2],decreasing=T),1:3]
ri.sort

ratio.final.phase.partial.pred.med <- pred.med
ratio.final.phase.partial.pred.up <- pred.up
ratio.final.phase.partial.pred.lo <- pred.lo
ratio.final.phase.partial.val.med <- val.med

# write outputs
dir.path <- paste(getwd(),"/out/",sep="")
for (o in 1:length(predictors.sub)) {
  var.name <- colnames(final.dat.imputed.orig.scale)[predictors.sub[o]]
  var.name.abbr <- toupper(sub("\\..*", "", var.name, ignore.case = T))
  var.title <- paste("ratio.final.phase.partial.",var.name.abbr,sep="")
  assign(var.title, 
         data.frame(val=ratio.final.phase.partial.val.med[,o],med=as.numeric(ratio.final.phase.partial.pred.med[,o]), 
                    up=as.numeric(ratio.final.phase.partial.pred.up[,o]), lo=as.numeric(ratio.final.phase.partial.pred.lo[,o])))
  get(var.title)
  write.table(get(var.title),paste(dir.path,paste(var.title,".csv",sep=""),sep=""), sep=",", row.names=F, col.names=T)
}

var.names <- rownames(ri.sort)
var.names.abbr <- toupper(sub("\\..*", "", var.names, ignore.case = T))
rel.infl.out <- data.frame(var=var.names.abbr, med=ri.sort$ri.med, up=ri.sort$ri.up, lo=ri.sort$ri.lo)
write.table(rel.infl.out,paste(dir.path,"ratio.final.phase.rel.infl.csv", sep=""), sep=",", row.names = F, col.names = T)



##################################
## general least-squares models ##
##################################

# get world map
wmap <- getMap(resolution="high")

# country centroids 
# calculate # centroids in a global equal-area CRS, then transform coordinates back to WGS84
wmap.sf <- sf::st_as_sf(wmap)
centroids.sf <- sf::st_transform(
  sf::st_centroid(sf::st_transform(wmap.sf, "ESRI:54034")),
  4326
)
centroid.coordinates <- sf::st_coordinates(centroids.sf)
centroids.df <- data.frame(
  lon=centroid.coordinates[, "X"],
  lat=centroid.coordinates[, "Y"],
  country=centroids.sf$NAME
)
head(centroids.df)

# response = ratio
final.dat.imputed.orig.scale$cntry.code
cc.lab <- data.frame(cont.cntry$country,cont.cntry$cntry.code)
colnames(cc.lab) <- c("country","cntry.code")
head(cc.lab)
final.dat.imputed.orig.scale2 <- merge(final.dat.imputed.orig.scale, cc.lab, by="cntry.code")
final.dat.imputed.orig.scale3 <- final.dat.imputed.orig.scale2[!duplicated(final.dat.imputed.orig.scale2$cntry.code), ]
head(final.dat.imputed.orig.scale3)
dim(final.dat.imputed.orig.scale3)

final.dat.imputed.orig.scale4 <- merge(final.dat.imputed.orig.scale3, centroids.df, by="country")
head(final.dat.imputed.orig.scale4)
tail(final.dat.imputed.orig.scale4)

x.sub <- which(colnames(final.dat.imputed.orig.scale4)=='lon')
y.sub <- which(colnames(final.dat.imputed.orig.scale4)=='lat')

final.dat.imputed.orig.scale4$x <- latlong2grid(final.dat.imputed.orig.scale4[,c(x.sub,y.sub)])$x # equidistant coordinates
final.dat.imputed.orig.scale4$y <- latlong2grid(final.dat.imputed.orig.scale4[,c(x.sub,y.sub)])$y # equidistant coordinates
plot(final.dat.imputed.orig.scale4$x,final.dat.imputed.orig.scale4$y,pch=19)

## determine best correlation structure
# ghsi, igs, agrL
colnames(final.dat.imputed.orig.scale4)

global.gls.formula <- ratio.l ~
  ghsi.sqr +
  igs.l +
  cpi.sqr +
  agrL.l +
  govexpedu

m1 <- gls(
  global.gls.formula,
  method = "ML",
  data = final.dat.imputed.orig.scale4
)

m2 <- gls(
  global.gls.formula,
  correlation = corExp(
    form = ~ lon + lat,
    nugget = TRUE
  ),
  method = "ML",
  data = final.dat.imputed.orig.scale4
)

m3 <- gls(
  global.gls.formula,
  correlation = corGaus(
    form = ~ lon + lat,
    nugget = TRUE
  ),
  method = "ML",
  data = final.dat.imputed.orig.scale4
)

m4 <- gls(
  global.gls.formula,
  correlation = corSpher(
    form = ~ lon + lat,
    nugget = TRUE
  ),
  method = "ML",
  data = final.dat.imputed.orig.scale4
)

m5 <- gls(
  global.gls.formula,
  correlation = corRatio(
    form = ~ lon + lat,
    nugget = TRUE
  ),
  method = "ML",
  data = final.dat.imputed.orig.scale4
)

vario1 <- Variogram(m1, form = ~ lon + lat, resType = "pearson")
plot(vario1, smooth = TRUE)

m2 <- gls(ratio.l ~ ghsi.sqr + igs.l + agrL.l, correlation = corExp(form = ~ lon + lat, nugget = T), data = final.dat.imputed.orig.scale4)
m3 <- gls(ratio.l ~ ghsi.sqr + igs.l + agrL.l, correlation = corGaus(form = ~ lon + lat, nugget = T), data = final.dat.imputed.orig.scale4)
m4 <- gls(ratio.l ~ ghsi.sqr + igs.l + agrL.l, correlation = corSpher(form = ~ lon + lat, nugget = T), data = final.dat.imputed.orig.scale4)
m5 <- gls(ratio.l ~ ghsi.sqr + igs.l + agrL.l, correlation = corRatio(form = ~ lon + lat, nugget = T), data = final.dat.imputed.orig.scale4)

mod.lab <- c("noCor","Exp","Gaus","Spher","Ratio")
AIC.vec <- c(AICc(m1), AICc(m2), AICc(m3), AICc(m4), AICc(m5))
dAIC.vec <- delta.IC(AIC.vec)
wAIC.vec <- weight.IC(dAIC.vec)
psR2.mcf <- c(nagelkerke(m1)$Pseudo.R.squared.for.model.vs.null[1], nagelkerke(m2)$Pseudo.R.squared.for.model.vs.null[1], nagelkerke(m3)$Pseudo.R.squared.for.model.vs.null[1], nagelkerke(m4)$Pseudo.R.squared.for.model.vs.null[1], nagelkerke(m5)$Pseudo.R.squared.for.model.vs.null[1])
psR2.cs <- c(nagelkerke(m1)$Pseudo.R.squared.for.model.vs.null[2], nagelkerke(m2)$Pseudo.R.squared.for.model.vs.null[2], nagelkerke(m3)$Pseudo.R.squared.for.model.vs.null[2], nagelkerke(m4)$Pseudo.R.squared.for.model.vs.null[2], nagelkerke(m5)$Pseudo.R.squared.for.model.vs.null[2])
psR2.cu <- c(nagelkerke(m1)$Pseudo.R.squared.for.model.vs.null[3], nagelkerke(m2)$Pseudo.R.squared.for.model.vs.null[3], nagelkerke(m3)$Pseudo.R.squared.for.model.vs.null[3], nagelkerke(m4)$Pseudo.R.squared.for.model.vs.null[3], nagelkerke(m5)$Pseudo.R.squared.for.model.vs.null[3])
results.out <- data.frame(mod.lab,AIC.vec,dAIC.vec,wAIC.vec,psR2.mcf,psR2.cs,psR2.cu)
colnames(results.out) <- c("mod","AICc","dAICc","wAICc","psR2mcf","psR2cs","psR2cu")
results.sort <- results.out[order(results.out[,4],decreasing=T),1:7]
results.sort

# percentage of variance explained by geographic coordinates
100*(1 - results.sort[3,5]/(results.sort[1,5])) # pseudo R2 - McFadden
100*(1 - results.sort[3,6]/(results.sort[1,6])) # pseudo R2 - Cox & Snell
100*(1 - results.sort[3,7]/(results.sort[1,7])) # pseudo R2 - Craig & Uhler

vario3 <- Variogram(m3, form = ~ lon + lat, resType = "pearson")
plot(vario3, smooth = TRUE)
vario3.nr <- Variogram(m3, form = ~ lon + lat, resType = "normalized")
plot(vario3.nr, smooth = TRUE)

# run with Gaussian spatial autocorrelation
# model set
vars <- c("ghsi.sqr","igs.l","agrL.l")    
vars3comb <- paste("ratio.l~", apply(combn(vars,3),2,paste,collapse='+'), sep="")
vars2comb <- paste("ratio.l~", apply(combn(vars,2),2,paste,collapse='+'), sep="")
vars1comb <- paste("ratio.l~", vars, sep="")
interconly <- "ratio.l~1"
mod.vec <- c(vars3comb,vars2comb,vars1comb,interconly)
mod.vec  

## define n.mod
n.mod <- length(mod.vec)

# model fitting and logLik output loop
Modnum <- length(mod.vec)
SaveCount <- BIC.vec <- AICc.vec <- LL.vec <- k.vec <- psR2mcf <- psR2cs <- psR2cu <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  #fit <- gls(as.formula(mod.vec[i]), correlation = corGaus(form = ~ lon + lat, nugget = T), method="ML", data = final.dat.imputed.orig.scale4)
  fit <- gls(as.formula(mod.vec[i]), method="ML", data = final.dat.imputed.orig.scale4)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- logLik(fit)
  k.vec[i] <-  (AIC(fit) - -2*as.numeric(logLik(fit)))/2
  AICc.vec[i] <- AICc(fit)
  BIC.vec[i] <- BIC(fit)
  psR2mcf[i] <- nagelkerke(fit)$Pseudo.R.squared.for.model.vs.null[1] # McFadden pseudo-R2
  psR2cs[i] <- nagelkerke(fit)$Pseudo.R.squared.for.model.vs.null[2] # Cox & Snell pseudo-R2
  psR2cu[i] <- nagelkerke(fit)$Pseudo.R.squared.for.model.vs.null[3] # Craig & Uhler pseudo-R2
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

# re-interpret pseudo-R2s as % deviance from maximum pseudo-R2
psR2mcf.rel <- 100 - 100*(max(psR2mcf) - psR2mcf)/max(psR2mcf)
psR2cs.rel <- 100 - 100*(max(psR2cs) - psR2cs)/max(psR2cs)
psR2cu.rel <- 100 - 100*(max(psR2cu) - psR2cu)/max(psR2cu)

sumtableGLS1 <- data.frame(mod.num,k.vec,round(LL.vec,3),round(AICc.vec,3),round(dAICc,3),round(wAICc,3),BIC.vec,round(dBIC,3),round(wBIC,3),round(psR2mcf,3),round(psR2mcf.rel,1),round(psR2cs,3),round(psR2cs.rel,1),round(psR2cu,3),round(psR2cu.rel,1))
colnames(sumtableGLS1) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","psR2mcf","mcfRel","psR2cs","csRel","psR2cu","cuRel")
row.names(sumtableGLS1) <- as.character(mod.vec)
summary.tableGLS1 <- sumtableGLS1[order(sumtableGLS1[,6],decreasing=T),] # order by wBIC
summary.tableGLS1

top.mod.subGLS1 <- summary.tableGLS1[1,1]
summary(mod.list[[top.mod.subGLS1]])
variotopGLS1 <- Variogram(mod.list[[top.mod.subGLS1]], form = ~ lon + lat, resType = "pearson")
plot(variotopGLS1, smooth = F, pch=19, grid=T)
variotop.nrGLS1 <- Variogram(mod.list[[top.mod.subGLS1]], form = ~ lon + lat, resType = "normalized")
plot(variotop.nrGLS1, smooth = F, pch=19, grid=T)




####################
####################
## response: pdam
####################
####################
head(final.datPdam.imputed.orig.scale)

######################
# resampled BRT loop #
######################
biter <- 1000
eq.sp.points <- 100

## variable selection
colnames(final.datPdam.imputed.orig.scale)
pred.gdp.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='gdp.sqr')
pred.cpi.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='cpi.sqr')
pred.ghsi.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='ghsi.sqr')
pred.igs.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='igs.l')
pred.agrL.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='agrL.l')
pred.VAag.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='VAag.l')
pred.govexpedu.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='govexpedu')
pred.stja.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='stja.l')
predictors.sub <- c(pred.gdp.sub,pred.cpi.sub,pred.ghsi.sub,pred.igs.sub,pred.agrL.sub,
                    pred.VAag.sub,pred.govexpedu.sub,pred.stja.sub)
predictors.sub
resp.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='pdam.l')
resp.sub

# create storage arrays
val.arr <- pred.arr <- array(data = 0, dim = c(eq.sp.points, length(predictors.sub), biter),
                             dimnames=list(paste("x",1:eq.sp.points,sep=""),
                                           attr(final.datPdam.imputed.orig.scale, "names")[predictors.sub], paste("b",1:biter,sep="")))

# create storage vectors
D2.vec <- CV.cor.vec <- CV.cor.se.vec <- 
  GDP.ri <- CPI.ri <- GHSI.ri <- IGS.ri <- AGRL.ri <- VAAG.ri <- GOVEXPEDU.ri <- STJA.ri <- rep(0,biter)

for (b in 1:biter) {
  # resample data among countries
  resamp.sub <- sort(sample(x = 1:dim(final.datPdam.imputed.orig.scale)[1], size = dim(final.datPdam.imputed.orig.scale)[1], replace=TRUE))
  dat.resamp <- final.datPdam.imputed.orig.scale[resamp.sub,]
  
  # boosted regression tree
  brt.fit <- adaptive.gbm.step(dat.resamp, gbm.x = attr(dat.resamp, "names")[predictors.sub],
                               gbm.y = attr(dat.resamp, "names")[resp.sub], family="gaussian", 
                               max.trees=100000, tolerance = 0.0001, learning.rate = 0.001, 
                               bag.fraction=0.75, tree.complexity = 2, silent=T, tolerance.method = "auto")
  summ.fit <- summary(brt.fit)
  
  # variable relative importance
  GDP.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[1]])])
  CPI.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[2]])])
  GHSI.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[3]])])
  IGS.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[4]])])
  AGRL.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[5]])])
  VAAG.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[6]])])
  GOVEXPEDU.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[7]])])
  STJA.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[8]])])
  
  D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / 
    brt.fit$cv.statistics$deviance.mean
  D2.vec[b] <- D2
  CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
  CV.cor.vec[b] <- CV.cor
  CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
  CV.cor.se.vec[b] <- CV.cor.se
  
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=length(predictors.sub))
  ## output average predictions
  for (p in 1:length(predictors.sub)) {
    RESP.val[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,2]
  }
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.fit$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.fit$var.names
  
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  print(b)
  
} # end b

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- 5
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] > (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  }
  print(k)
}

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)

val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median)

## plot
# order: GDP,CPI,GHSI,IGS,AGRL,VAAG,GOVEXPEDU,STJA
par(mfrow=c(2,4)) 

plot(val.med[,1],pred.med[,1],type="l",ylim=c(min(pred.lo[,1]),max(pred.up[,1])), lwd=2,
     ylab="(←lower) proportion damage (higher→)", xlab="(←poorer) GDPpc (richer→)" )
lines(val.med[,1], pred.lo[,1], type="l", lty=2, col="red")
lines(val.med[,1], pred.up[,1], type="l", lty=2, col="red")

plot(val.med[,2],pred.med[,2],type="l",ylim=c(min(pred.lo[,2]),max(pred.up[,2])), lwd=2,
     ylab="(←lower) proportion damage (higher→)", xlab="(←higher) corruption (lower→)" )
lines(val.med[,2], pred.lo[,2], type="l", lty=2, col="red")
lines(val.med[,2], pred.up[,2], type="l", lty=2, col="red")

plot(val.med[,3],pred.med[,3],type="l",ylim=c(min(pred.lo[,3]),max(pred.up[,3])), lwd=2,
     ylab="(←lower) proportion damage (higher→)", xlab="(←lower) global health security (higher→)" )
lines(val.med[,3], pred.lo[,3], type="l", lty=3, col="red")
lines(val.med[,3], pred.up[,3], type="l", lty=3, col="red")

plot(val.med[,4],pred.med[,4],type="l",ylim=c(min(pred.lo[,4]),max(pred.up[,4])), lwd=2,
     ylab="(←lower) proportion damage (higher→)", xlab="(←less) goods & services imports (more→)" )
lines(val.med[,4], pred.lo[,4], type="l", lty=3, col="red")
lines(val.med[,4], pred.up[,4], type="l", lty=3, col="red")

plot(val.med[,5],pred.med[,5],type="l",ylim=c(min(pred.lo[,5]),max(pred.up[,5])), lwd=2,
     ylab="(←lower) proportion damage (higher→)", xlab="(←less) agricultural land (more→)" )
lines(val.med[,5], pred.lo[,5], type="l", lty=2, col="red")
lines(val.med[,5], pred.up[,5], type="l", lty=2, col="red")

plot(val.med[,6],pred.med[,6],type="l",ylim=c(min(pred.lo[,6]),max(pred.up[,6])), lwd=2,
     ylab="(←lower) proportion damage (higher→)", xlab="(←less) value added primary production (more→)" )
lines(val.med[,6], pred.lo[,6], type="l", lty=2, col="red")
lines(val.med[,6], pred.up[,6], type="l", lty=2, col="red")

plot(val.med[,7],pred.med[,7],type="l",ylim=c(min(pred.lo[,7]),max(pred.up[,7])), lwd=2,
     ylab="(←lower) proportion damage (higher→)", xlab="(←less) education investment (more→)" )
lines(val.med[,7], pred.lo[,7], type="l", lty=2, col="red")
lines(val.med[,7], pred.up[,7], type="l", lty=2, col="red")

plot(val.med[,8],pred.med[,8],type="l",ylim=c(min(pred.lo[,8]),max(pred.up[,8])), lwd=2,
     ylab="(←lower) proportion damage (higher→)", xlab="(←fewer) scientific & tech journal articles (more→)" )
lines(val.med[,8], pred.lo[,8], type="l", lty=2, col="red")
lines(val.med[,8], pred.up[,8], type="l", lty=2, col="red")

par(mfrow=c(1,1)) 

# kappa method for output vectors
D2.update <- D2.vec[1:biter]
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# order: GDP,CPI,GHSI,IGS,AGRL,VAAG,GOVEXPEDU,STJA
GDP.ri.update <- GDP.ri[1:biter]
CPI.ri.update <- CPI.ri[1:biter]
GHSI.ri.update <- GHSI.ri[1:biter]
IGS.ri.update <- IGS.ri[1:biter]
AGRL.ri.update <- AGRL.ri[1:biter]
VAAG.ri.update <- VAAG.ri[1:biter]
GOVEXPEDU.ri.update <- GOVEXPEDU.ri[1:biter]
STJA.ri.update <- STJA.ri[1:biter]

for (k in 1:kappa.n) {
  D2.mean <- mean(D2.update, na.rm=T); D2.sd <- sd(D2.update, na.rm=T)
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  GDP.mean <- mean(GDP.ri.update, na.rm=T); GDP.sd <- sd(GDP.ri.update, na.rm=T)
  CPI.mean <- mean(CPI.ri.update, na.rm=T); CPI.sd <- sd(CPI.ri.update, na.rm=T)
  GHSI.mean <- mean(GHSI.ri.update, na.rm=T); GHSI.sd <- sd(GHSI.ri.update, na.rm=T)
  IGS.mean <- mean(IGS.ri.update, na.rm=T); IGS.sd <- sd(IGS.ri.update, na.rm=T)
  AGRL.mean <- mean(AGRL.ri.update, na.rm=T); AGRL.sd <- sd(AGRL.ri.update, na.rm=T)
  VAAG.mean <- mean(VAAG.ri.update, na.rm=T); VAAG.sd <- sd(VAAG.ri.update, na.rm=T)
  GOVEXPEDU.mean <- mean(GOVEXPEDU.ri.update, na.rm=T); GOVEXPEDU.sd <- sd(GOVEXPEDU.ri.update, na.rm=T)
  STJA.mean <- mean(STJA.ri.update, na.rm=T); STJA.sd <- sd(STJA.ri.update, na.rm=T)
  
  for (u in 1:biter) {
    D2.update[u] <- ifelse((D2.update[u] < (D2.mean-kappa*D2.sd) | D2.update[u] > (D2.mean+kappa*D2.sd)), NA, D2.update[u])
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] > (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] > (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    
    # order: GDP,CPI,GHSI,IGS,AGRL,VAAG,GOVEXPEDU,STJA
    GDP.ri.update[u] <- ifelse((GDP.ri.update[u] < (GDP.mean-kappa*GDP.sd) | GDP.ri.update[u] > (GDP.mean+kappa*GDP.sd)), NA, GDP.ri.update[u])
    CPI.ri.update[u] <- ifelse((CPI.ri.update[u] < (CPI.mean-kappa*CPI.sd) | CPI.ri.update[u] > (CPI.mean+kappa*CPI.sd)), NA, CPI.ri.update[u])
    GHSI.ri.update[u] <- ifelse((GHSI.ri.update[u] < (GHSI.mean-kappa*GHSI.sd) | GHSI.ri.update[u] > (GHSI.mean+kappa*GHSI.sd)), NA, GHSI.ri.update[u])
    IGS.ri.update[u] <- ifelse((IGS.ri.update[u] < (IGS.mean-kappa*IGS.sd) | IGS.ri.update[u] > (IGS.mean+kappa*IGS.sd)), NA, IGS.ri.update[u])
    AGRL.ri.update[u] <- ifelse((AGRL.ri.update[u] < (AGRL.mean-kappa*AGRL.sd) | AGRL.ri.update[u] > (AGRL.mean+kappa*AGRL.sd)), NA, AGRL.ri.update[u])
    VAAG.ri.update[u] <- ifelse((VAAG.ri.update[u] < (VAAG.mean-kappa*VAAG.sd) | VAAG.ri.update[u] > (VAAG.mean+kappa*VAAG.sd)), NA, VAAG.ri.update[u])
    GOVEXPEDU.ri.update[u] <- ifelse((GOVEXPEDU.ri.update[u] < (GOVEXPEDU.mean-kappa*GOVEXPEDU.sd) | GOVEXPEDU.ri.update[u] > (GOVEXPEDU.mean+kappa*GOVEXPEDU.sd)), NA, GOVEXPEDU.ri.update[u])
    STJA.ri.update[u] <- ifelse((STJA.ri.update[u] < (STJA.mean-kappa*STJA.sd) | STJA.ri.update[u] > (STJA.mean+kappa*STJA.sd)), NA, STJA.ri.update[u])
  }
  print(k)
}

D2.med <- median(D2.update, na.rm=TRUE)
D2.lo <- quantile(D2.update, probs=0.025, na.rm=TRUE)
D2.up <- quantile(D2.update, probs=0.975, na.rm=TRUE)
print(c(D2.lo,D2.med,D2.up))

CV.cor.med <- median(CV.cor.update, na.rm=TRUE)
CV.cor.lo <- quantile(CV.cor.update, probs=0.025, na.rm=TRUE)
CV.cor.up <- quantile(CV.cor.update, probs=0.975, na.rm=TRUE)
print(c(CV.cor.lo,CV.cor.med,CV.cor.up))

# order: GDP,CPI,GHSI,IGS,AGRL,VAAG,GOVEXPEDU,STJA
GDP.ri.lo <- quantile(GDP.ri.update, probs=0.025, na.rm=TRUE)
GDP.ri.med <- median(GDP.ri.update, na.rm=TRUE)
GDP.ri.up <- quantile(GDP.ri.update, probs=0.975, na.rm=TRUE)

CPI.ri.lo <- quantile(CPI.ri.update, probs=0.025, na.rm=TRUE)
CPI.ri.med <- median(CPI.ri.update, na.rm=TRUE)
CPI.ri.up <- quantile(CPI.ri.update, probs=0.975, na.rm=TRUE)

GHSI.ri.lo <- quantile(GHSI.ri.update, probs=0.025, na.rm=TRUE)
GHSI.ri.med <- median(GHSI.ri.update, na.rm=TRUE)
GHSI.ri.up <- quantile(GHSI.ri.update, probs=0.975, na.rm=TRUE)

IGS.ri.lo <- quantile(IGS.ri.update, probs=0.025, na.rm=TRUE)
IGS.ri.med <- median(IGS.ri.update, na.rm=TRUE)
IGS.ri.up <- quantile(IGS.ri.update, probs=0.975, na.rm=TRUE)

AGRL.ri.lo <- quantile(AGRL.ri.update, probs=0.025, na.rm=TRUE)
AGRL.ri.med <- median(AGRL.ri.update, na.rm=TRUE)
AGRL.ri.up <- quantile(AGRL.ri.update, probs=0.975, na.rm=TRUE)

VAAG.ri.lo <- quantile(VAAG.ri.update, probs=0.025, na.rm=TRUE)
VAAG.ri.med <- median(VAAG.ri.update, na.rm=TRUE)
VAAG.ri.up <- quantile(VAAG.ri.update, probs=0.975, na.rm=TRUE)

GOVEXPEDU.ri.lo <- quantile(GOVEXPEDU.ri.update, probs=0.025, na.rm=TRUE)
GOVEXPEDU.ri.med <- median(GOVEXPEDU.ri.update, na.rm=TRUE)
GOVEXPEDU.ri.up <- quantile(GOVEXPEDU.ri.update, probs=0.975, na.rm=TRUE)

STJA.ri.lo <- quantile(STJA.ri.update, probs=0.025, na.rm=TRUE)
STJA.ri.med <- median(STJA.ri.update, na.rm=TRUE)
STJA.ri.up <- quantile(STJA.ri.update, probs=0.975, na.rm=TRUE)

# order: GDP,CPI,GHSI,IGS,AGRL,VAAG,GOVEXPEDU,STJA
ri.lo <- c(GDP.ri.lo,CPI.ri.lo,GHSI.ri.lo,IGS.ri.lo,AGRL.ri.lo,VAAG.ri.lo,GOVEXPEDU.ri.lo,STJA.ri.lo)
ri.med <- c(GDP.ri.med,CPI.ri.med,GHSI.ri.med,IGS.ri.med,AGRL.ri.med,VAAG.ri.med,GOVEXPEDU.ri.med,STJA.ri.med)
ri.up <- c(GDP.ri.up,CPI.ri.up,GHSI.ri.up,IGS.ri.up,AGRL.ri.up,VAAG.ri.up,GOVEXPEDU.ri.up,STJA.ri.up)

ri.out <- as.data.frame(cbind(ri.lo,ri.med,ri.up))
colnames(ri.out) <- c("ri.lo","ri.med","ri.up")
rownames(ri.out) <- attr(final.datPdam.imputed.orig.scale, "names")[predictors.sub]
ri.sort <- ri.out[order(ri.out[,2],decreasing=T),1:3]
ri.sort

pdam.partial.pred.med <- pred.med
pdam.partial.pred.up <- pred.up
pdam.partial.pred.lo <- pred.lo
pdam.partial.val.med <- val.med

# write outputs
dir.path <- paste(getwd(),"/out/",sep="")
for (o in 1:length(predictors.sub)) {
  var.name <- colnames(final.dat.imputed.orig.scale)[predictors.sub[o]]
  var.name.abbr <- toupper(sub("\\..*", "", var.name, ignore.case = T))
  var.title <- paste("pdam.partial.",var.name.abbr,sep="")
  assign(var.title, 
         data.frame(val=pdam.partial.val.med[,o],med=as.numeric(pdam.partial.pred.med[,o]), 
                    up=as.numeric(pdam.partial.pred.up[,o]), lo=as.numeric(pdam.partial.pred.lo[,o])))
  get(var.title)
  write.table(get(var.title),paste(dir.path,paste(var.title,".csv",sep=""),sep=""), sep=",", row.names=F, col.names=T)
}

var.names <- rownames(ri.sort)
var.names.abbr <- toupper(sub("\\..*", "", var.names, ignore.case = T))
rel.infl.out <- data.frame(var=var.names.abbr, med=ri.sort$ri.med, up=ri.sort$ri.up, lo=ri.sort$ri.lo)
write.table(rel.infl.out,paste(dir.path,"pdam.rel.infl.csv", sep=""), sep=",", row.names = F, col.names = T)



## PHASE 1 - wealth/capacity
# gdp, cpi, ghsi, igs
colnames(final.dat.imputed.orig.scale)
pred.gdp.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='gdp.sqr')
pred.cpi.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='cpi.sqr')
pred.ghsi.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='ghsi.sqr')
pred.igs.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='igs.l')
predictors.sub <- c(pred.gdp.sub,pred.cpi.sub,pred.ghsi.sub,pred.igs.sub)
resp.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='pdam.l')

brt.fit4 <- adaptive.gbm.step(final.datPdam.imputed.orig.scale, gbm.x = attr(final.datPdam.imputed.orig.scale, "names")[predictors.sub],
                              gbm.y = attr(final.datPdam.imputed.orig.scale, "names")[resp.sub], family="gaussian", max.trees=100000, tolerance = 0.00003, learning.rate = 0.00001, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit4)
D2 <- 100 * (brt.fit4$cv.statistics$deviance.mean - brt.fit4$self.statistics$mean.resid) / brt.fit4$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit4)
gbm.plot.fits(brt.fit4)

brt4.CV.cor <- 100 * brt.fit4$cv.statistics$correlation.mean
brt4.CV.cor
brt4.CV.cor.se <- 100 * brt.fit4$cv.statistics$correlation.se
brt4.CV.cor.se
print(c(brt4.CV.cor, brt4.CV.cor.se))


# GLMM
# recode continent to greater region to increase per-category sample sizes
table(final.datPdam.imputed.orig.scale$cont)
final.datPdam.imputed.orig.scale$reg2 <- ifelse(final.datPdam.imputed.orig.scale$cont == "NAM" | final.datPdam.imputed.orig.scale$cont == "CAR", "NAMCAR", final.datPdam.imputed.orig.scale$cont)
final.datPdam.imputed.orig.scale$reg2 <- ifelse(final.datPdam.imputed.orig.scale$cont == "EUR" | final.datPdam.imputed.orig.scale$cont == "ME", "EURME", final.datPdam.imputed.orig.scale$reg2)
final.datPdam.imputed.orig.scale$reg2 <- ifelse(final.datPdam.imputed.orig.scale$cont == "ASIA" | final.datPdam.imputed.orig.scale$cont == "OC", "ASIAOC", final.datPdam.imputed.orig.scale$reg2)
table(final.datPdam.imputed.orig.scale$reg2)

# model set
colnames(final.datPdam.imputed.orig.scale)
vars <- c("gdp.sqr","cpi.sqr","ghsi.sqr","igs.l")    
vars4comb <- paste("pdam.l~", apply(combn(vars,4),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars3comb <- paste("pdam.l~", apply(combn(vars,3),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars2comb <- paste("pdam.l~", apply(combn(vars,2),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars1comb <- paste("pdam.l~", vars, "+(1|reg2)", sep="")
interconly <- "pdam.l~1+(1|reg2)"
mod.vec <- c(vars4comb,vars3comb,vars2comb,vars1comb,interconly)
mod.vec  

length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.datPdam.imputed.orig.scale, na.action=na.omit, REML=FALSE, verbose=F)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- AICc(fit)
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r2_nakagawa(fit)$R2_marginal # marginal R-squared
  Rc[i] <- 100*r2_nakagawa(fit)$R2_conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable1 <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable1) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable1) <- as.character(mod.vec)
summary.table1 <- sumtable1[order(sumtable1[,5],decreasing=F),]
summary.table1

fitsat <- lmer(as.formula(mod.vec[1]), data=final.datPdam.imputed.orig.scale, na.action=na.omit, REML=FALSE)
check_model(fitsat, detrend=T)
plot_model(fitsat, show.values=T, vline.color = "purple")
vif(fitsat)

## gdp collinear; remove
# model set
colnames(final.datPdam.imputed.orig.scale)
vars <- c("cpi.sqr","ghsi.sqr","igs.l")    
vars3comb <- paste("pdam.l~", apply(combn(vars,3),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars2comb <- paste("pdam.l~", apply(combn(vars,2),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars1comb <- paste("pdam.l~", vars, "+(1|reg2)", sep="")
interconly <- "pdam.l~1+(1|reg2)"
mod.vec <- c(vars3comb,vars2comb,vars1comb,interconly)
mod.vec  

length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.datPdam.imputed.orig.scale, na.action=na.omit, REML=FALSE, verbose=F)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- AICc(fit)
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r2_nakagawa(fit)$R2_marginal # marginal R-squared
  Rc[i] <- 100*r2_nakagawa(fit)$R2_conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable1 <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable1) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable1) <- as.character(mod.vec)
summary.table1 <- sumtable1[order(sumtable1[,5],decreasing=F),]
summary.table1

fitsat <- lmer(as.formula(mod.vec[1]), data=final.datPdam.imputed.orig.scale, na.action=na.omit, REML=FALSE)
check_model(fitsat, detrend=T)
plot_model(fitsat, show.values=T, vline.color = "purple")
vif(fitsat)


## PHASE 2 - primary production
# agrL, VAag
colnames(final.datPdam.imputed.orig.scale)
pred.agrl.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='agrL.l')
pred.vaag.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='VAag.l')
predictors.sub <- c(pred.agrl.sub,pred.vaag.sub)
resp.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='pdam.l')

brt.fit1 <- adaptive.gbm.step(final.datPdam.imputed.orig.scale, gbm.x = attr(final.datPdam.imputed.orig.scale, "names")[predictors.sub],
                              gbm.y = attr(final.datPdam.imputed.orig.scale, "names")[resp.sub], family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.00001, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit1)
D2 <- 100 * (brt.fit1$cv.statistics$deviance.mean - brt.fit1$self.statistics$mean.resid) / brt.fit1$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit1)
gbm.plot.fits(brt.fit1)

brt1.CV.cor <- 100 * brt.fit1$cv.statistics$correlation.mean
brt1.CV.cor
brt1.CV.cor.se <- 100 * brt.fit1$cv.statistics$correlation.se
brt1.CV.cor.se
print(c(brt1.CV.cor, brt1.CV.cor.se))

# model set
vars <- c("agrL.l","VAag.l")    
vars2comb <- paste("pdam.l~", apply(combn(vars,2),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars1comb <- paste("pdam.l~", vars, "+(1|reg2)", sep="")
interconly <- "pdam.l~1+(1|reg2)"
mod.vec <- c(vars2comb,vars1comb,interconly)
mod.vec  

length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.datPdam.imputed.orig.scale, na.action=na.omit, REML=FALSE)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- AICc(fit)
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r2_nakagawa(fit)$R2_marginal # marginal R-squared
  Rc[i] <- 100*r2_nakagawa(fit)$R2_conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable2 <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable2) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable2) <- as.character(mod.vec)
summary.table2 <- sumtable2[order(sumtable2[,5],decreasing=F),]
summary.table2

fitsat <- lmer(as.formula(mod.vec[1]), data=final.datPdam.imputed.orig.scale, na.action=na.omit, REML=FALSE)
check_model(fitsat, detrend=T)
plot_model(fitsat, show.values=T, vline.color = "purple")


## PHASE 3 - education/research
# govexpedu, stja
colnames(final.datPdam.imputed.orig.scale)
pred.govexpedu.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='govexpedu')
pred.stja.l.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='stja.l')
predictors.sub <- c(pred.govexpedu.sub,pred.stja.l.sub)
resp.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='pdam.l')

brt.fit1 <- adaptive.gbm.step(final.datPdam.imputed.orig.scale, gbm.x = attr(final.datPdam.imputed.orig.scale, "names")[predictors.sub],
                              gbm.y = attr(final.datPdam.imputed.orig.scale, "names")[resp.sub], family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.00001, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit1)
D2 <- 100 * (brt.fit1$cv.statistics$deviance.mean - brt.fit1$self.statistics$mean.resid) / brt.fit1$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit1)
gbm.plot.fits(brt.fit1)

brt1.CV.cor <- 100 * brt.fit1$cv.statistics$correlation.mean
brt1.CV.cor
brt1.CV.cor.se <- 100 * brt.fit1$cv.statistics$correlation.se
brt1.CV.cor.se
print(c(brt1.CV.cor, brt1.CV.cor.se))


# GLMM
# model set
colnames(final.datPdam.imputed.orig.scale)
vars <- c("govexpedu","stja.l")    
vars2comb <- paste("pdam.l~", apply(combn(vars,2),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars1comb <- paste("pdam.l~", vars, "+(1|reg2)", sep="")
interconly <- "pdam.l~1+(1|reg2)"
mod.vec <- c(vars2comb,vars1comb,interconly)
mod.vec  

length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.datPdam.imputed.orig.scale, na.action=na.omit, REML=FALSE)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- AIC(fit)
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r2_nakagawa(fit)$R2_marginal # marginal R-squared
  Rc[i] <- 100*r2_nakagawa(fit)$R2_conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable3 <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable3) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable3) <- as.character(mod.vec)
summary.table3 <- sumtable3[order(sumtable3[,5],decreasing=F),]
summary.table3

fitsat <- lmer(as.formula(mod.vec[1]), data=final.datPdam.imputed.orig.scale, na.action=na.omit, REML=FALSE)
check_model(fitsat, detrend=T)
plot_model(fitsat, show.values=T, vline.color = "purple")



## COMBINE PHASES
# igs, ghsi, igs, VAag, stja, govexpedu
# model set
colnames(final.datPdam.imputed.orig.scale)
vars <- c("cpi.sqr","ghsi.sqr","igs.l","VAag.l","stja.l","govexpedu")    
vars6comb <- paste("pdam.l~", apply(combn(vars,6),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars5comb <- paste("pdam.l~", apply(combn(vars,5),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars4comb <- paste("pdam.l~", apply(combn(vars,4),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars3comb <- paste("pdam.l~", apply(combn(vars,3),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars2comb <- paste("pdam.l~", apply(combn(vars,2),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars1comb <- paste("pdam.l~", vars, "+(1|reg2)", sep="")
interconly <- "pdam.l~1+(1|reg2)"
mod.vec <- c(vars6comb,vars5comb,vars4comb,vars3comb,vars2comb,vars1comb,interconly)
mod.vec  

length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.datPdam.imputed.orig.scale, na.action=na.omit, REML=FALSE)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- AICc(fit)
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r2_nakagawa(fit)$R2_marginal # marginal R-squared
  Rc[i] <- 100*r2_nakagawa(fit)$R2_conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtableF <- data.frame(mod.num,k.vec,round(LL.vec,3),round(AICc.vec,3),round(dAICc,3),round(wAICc,3),
                        round(BIC.vec,3),round(dBIC,3),round(wBIC,4),round(Rm,1),round(Rc,1))
colnames(sumtableF) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtableF) <- as.character(mod.vec)
summary.tableF <- sumtableF[order(sumtableF[,5],decreasing=F),]
summary.tableF

fitsat <- lmer(as.formula(mod.vec[1]), data=final.datPdam.imputed.orig.scale, na.action=na.omit, REML=FALSE)
vif(fitsat)
check_model(fitsat, detrend=F)
plot_model(fitsat, type="est", sort.est=T)


## BRT
# cpi, ghsi, igs, VAag, stja, govexpedu
colnames(final.datPdam.imputed.orig.scale)
pred.cpi.sub <- which(colnames((final.datPdam.imputed.orig.scale))=='cpi.sqr')
pred.ghsi.sub <- which(colnames((final.datPdam.imputed.orig.scale))=='ghsi.sqr')
pred.igs.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='igs.l')
pred.VAag.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='VAag.l')
pred.stja.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='stja.l')
pred.govexpedu.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='govexpedu')
predictors.sub <- c(pred.cpi.sub,pred.ghsi.sub,pred.igs.sub,pred.VAag.sub,pred.stja.sub,pred.govexpedu.sub)
resp.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='pdam.l')

brt.fitF <- adaptive.gbm.step(final.datPdam.imputed.orig.scale, gbm.x = attr(final.datPdam.imputed.orig.scale, "names")[predictors.sub],
                              gbm.y = attr(final.datPdam.imputed.orig.scale, "names")[resp.sub], family="gaussian", 
                              max.trees=100000, tolerance = 0.00001, learning.rate = 0.00003, 
                              bag.fraction=0.75, tree.complexity = 2)
summary(brt.fitF)
D2 <- 100 * (brt.fitF$cv.statistics$deviance.mean - brt.fitF$self.statistics$mean.resid) / brt.fitF$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fitF)
gbm.plot.fits(brt.fitF)

brtF.CV.cor <- 100 * brt.fitF$cv.statistics$correlation.mean
brtF.CV.cor
brtF.CV.cor.se <- 100 * brt.fitF$cv.statistics$correlation.se
brtF.CV.cor.se
print(c(brtF.CV.cor, brtF.CV.cor.se))


# resampled BRT loop
# cpi, ghsi, igs, VAag, govexpedu
biter <- 200
eq.sp.points <- 100

colnames(final.datPdam.imputed.orig.scale)
pred.cpi.sub <- which(colnames((final.datPdam.imputed.orig.scale))=='cpi.sqr')
pred.ghsi.sub <- which(colnames((final.datPdam.imputed.orig.scale))=='ghsi.sqr')
pred.igs.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='igs.l')
pred.VAag.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='VAag.l')
pred.govexpedu.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='govexpedu')
predictors.sub <- c(pred.cpi.sub,pred.ghsi.sub,pred.igs.sub,pred.VAag.sub,pred.govexpedu.sub)
resp.sub <- which(colnames(final.datPdam.imputed.orig.scale)=='pdam.l')

# create storage arrays
val.arr <- pred.arr <- array(data = NA, dim = c(eq.sp.points, length(predictors.sub), biter),
                             dimnames=list(paste("x",1:eq.sp.points,sep=""), 
                                           attr(final.datPdam.imputed.orig.scale, "names")[predictors.sub], 
                                           paste("b",1:biter,sep="")))

# create storage vectors
D2.vec <- CV.cor.vec <- CV.cor.se.vec <- CPI.ri <- GHSI.ri <- IGS.ri <- VAAG.ri <- STJA.ri <- GOVEXPEDU.ri <- rep(NA,biter)

for (b in 1:biter) {
  # resample data among countries
  resamp.sub <- sort(sample(x = 1:dim(final.datPdam.imputed.orig.scale)[1], size = dim(final.datPdam.imputed.orig.scale)[1], replace=TRUE))
  dat.resamp <- final.datPdam.imputed.orig.scale[resamp.sub,]
  
  # boosted regression tree
  brt.fit <- adaptive.gbm.step(dat.resamp, gbm.x = attr(dat.resamp, "names")[predictors.sub],
                               gbm.y = attr(dat.resamp, "names")[resp.sub], family="gaussian", 
                               max.trees=100000, tolerance = 0.0001, learning.rate = 0.0001, 
                               bag.fraction=0.75, tree.complexity = 2, silent=T, tolerance.method = "auto")
  summ.fit <- summary(brt.fit)
  length(summ.fit[[1]])
  
  # variable relative importance
  CPI.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[1]])])
  GHSI.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[2]])])
  IGS.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[3]])])
  VAAG.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[4]])])
  GOVEXPEDU.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[5]])])

  
  D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
  D2.vec[b] <- D2
  CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
  CV.cor.vec[b] <- CV.cor
  CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
  CV.cor.se.vec[b] <- CV.cor.se
  
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=length(predictors.sub))
  ## output average predictions
  for (p in 1:length(predictors.sub)) {
    RESP.val[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,2]
  }
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.fit$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.fit$var.names
  
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  print(b)
  
} # end b

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- 5
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  }
  print(k)
}

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)

val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

par(mfrow=c(2,3)) 
# cpi, ghsi, igs, VAag, stja
plot(val.med[,1],pred.med[,1],type="l",ylim=c(min(pred.lo[,1]),max(pred.up[,1])),
     lwd=2, ylab="(←lower) proportion damage (higher→)",  xlab="(←higher) corruption (lower→)" )
lines(val.med[,1], pred.lo[,1], type="l", lty=2, col="red")
lines(val.med[,1], pred.up[,1], type="l", lty=2, col="red")

plot(val.med[,2],pred.med[,2],type="l",ylim=c(min(pred.lo[,2]),max(pred.up[,2])),
     lwd=2, ylab="(←lower) proportion damage (higher→)",  xlab="(←lower) global health security index (higher→)" )
lines(val.med[,2], pred.lo[,2], type="l", lty=2, col="red")
lines(val.med[,2], pred.up[,2], type="l", lty=2, col="red")

plot(val.med[,3],pred.med[,3],type="l",ylim=c(min(pred.lo[,3]),max(pred.up[,3])),
     lwd=2, ylab="(←lower) proportion damage (higher→)",  xlab="(←lower) per capita imports goods & services (higher→)" )
lines(val.med[,3], pred.lo[,3], type="l", lty=2, col="red")
lines(val.med[,3], pred.up[,3], type="l", lty=2, col="red")

plot(val.med[,4],pred.med[,4],type="l",ylim=c(min(pred.lo[,4]),max(pred.up[,4])),
     lwd=3, ylab="(←lower) proportion damage (higher→)", xlab="(←less) value-added agriculture (more→)")
lines(val.med[,4], pred.lo[,4], type="l", lty=2, col="red")
lines(val.med[,4], pred.up[,4], type="l", lty=2, col="red")

plot(val.med[,5],pred.med[,5],type="l",ylim=c(min(pred.lo[,5]),max(pred.up[,5])),
     lwd=3, ylab="(←lower) proportion damage (higher→)", xlab="(←lower) gov education expenditure (more→)")
lines(val.med[,5], pred.lo[,5], type="l", lty=2, col="red")
lines(val.med[,5], pred.up[,5], type="l", lty=2, col="red")

par(mfrow=c(1,1)) 

# kappa method for output vectors
D2.update <- D2.vec[1:biter]
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# igs, agrL, govexpedu
CPI.ri.update <- CPI.ri[1:biter]
GHSI.ri.update <- GHSI.ri[1:biter]
IGS.ri.update <- IGS.ri[1:biter]
VAAG.ri.update <- VAAG.ri[1:biter]
GOVEXPEDU.ri.update <- GOVEXPEDU.ri[1:biter]

for (k in 1:kappa.n) {
  D2.mean <- mean(D2.update, na.rm=T); D2.sd <- sd(D2.update, na.rm=T)
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  CPI.mean <- mean(CPI.ri.update, na.rm=T); CPI.sd <- sd(CPI.ri.update, na.rm=T)
  GHSI.mean <- mean(GHSI.ri.update, na.rm=T); GHSI.sd <- sd(GHSI.ri.update, na.rm=T)
  IGS.mean <- mean(IGS.ri.update, na.rm=T); IGS.sd <- sd(IGS.ri.update, na.rm=T)
  VAAG.mean <- mean(VAAG.ri.update, na.rm=T); VAAG.sd <- sd(VAAG.ri.update, na.rm=T)
  GOVEXPEDU.mean <- mean(GOVEXPEDU.ri.update, na.rm=T); GOVEXPEDU.sd <- sd(GOVEXPEDU.ri.update, na.rm=T)
  
  for (u in 1:biter) {
    D2.update[u] <- ifelse((D2.update[u] < (D2.mean-kappa*D2.sd) | D2.update[u] > (D2.mean+kappa*D2.sd)), NA, D2.update[u])
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] > (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] > (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    
    CPI.ri.update[u] <- ifelse((CPI.ri.update[u] < (CPI.mean-kappa*CPI.sd) | CPI.ri.update[u] > (CPI.mean+kappa*CPI.sd)), NA, CPI.ri.update[u])
    GHSI.ri.update[u] <- ifelse((GHSI.ri.update[u] < (GHSI.mean-kappa*GHSI.sd) | GHSI.ri.update[u] > (GHSI.mean+kappa*GHSI.sd)), NA, GHSI.ri.update[u])
    IGS.ri.update[u] <- ifelse((IGS.ri.update[u] < (IGS.mean-kappa*IGS.sd) | IGS.ri.update[u] > (IGS.mean+kappa*IGS.sd)), NA, IGS.ri.update[u])
    VAAG.ri.update[u] <- ifelse((VAAG.ri.update[u] < (VAAG.mean-kappa*VAAG.sd) | VAAG.ri.update[u] > (VAAG.mean+kappa*VAAG.sd)), NA, VAAG.ri.update[u])
    GOVEXPEDU.ri.update[u] <- ifelse((GOVEXPEDU.ri.update[u] < (GOVEXPEDU.mean-kappa*GOVEXPEDU.sd) | GOVEXPEDU.ri.update[u] > (GOVEXPEDU.mean+kappa*GOVEXPEDU.sd)), NA, GOVEXPEDU.ri.update[u])
  }
  
  print(k)
}

D2.med <- median(D2.update, na.rm=TRUE)
D2.lo <- quantile(D2.update, probs=0.025, na.rm=TRUE)
D2.up <- quantile(D2.update, probs=0.975, na.rm=TRUE)
print(c(D2.lo,D2.med,D2.up))

CV.cor.med <- median(CV.cor.update, na.rm=TRUE)
CV.cor.lo <- quantile(CV.cor.update, probs=0.025, na.rm=TRUE)
CV.cor.up <- quantile(CV.cor.update, probs=0.975, na.rm=TRUE)
print(c(CV.cor.lo,CV.cor.med,CV.cor.up))

CPI.ri.lo <- quantile(CPI.ri.update, probs=0.025, na.rm=TRUE)
CPI.ri.med <- median(CPI.ri.update, na.rm=TRUE)
CPI.ri.up <- quantile(CPI.ri.update, probs=0.975, na.rm=TRUE)

GHSI.ri.lo <- quantile(GHSI.ri.update, probs=0.025, na.rm=TRUE)
GHSI.ri.med <- median(GHSI.ri.update, na.rm=TRUE)
GHSI.ri.up <- quantile(GHSI.ri.update, probs=0.975, na.rm=TRUE)

IGS.ri.lo <- quantile(IGS.ri.update, probs=0.025, na.rm=TRUE)
IGS.ri.med <- median(IGS.ri.update, na.rm=TRUE)
IGS.ri.up <- quantile(IGS.ri.update, probs=0.975, na.rm=TRUE)

VAAG.ri.lo <- quantile(VAAG.ri.update, probs=0.025, na.rm=TRUE)
VAAG.ri.med <- median(VAAG.ri.update, na.rm=TRUE)
VAAG.ri.up <- quantile(VAAG.ri.update, probs=0.975, na.rm=TRUE)

GOVEXPEDU.ri.lo <- quantile(GOVEXPEDU.ri.update, probs=0.025, na.rm=TRUE)
GOVEXPEDU.ri.med <- median(GOVEXPEDU.ri.update, na.rm=TRUE)
GOVEXPEDU.ri.up <- quantile(GOVEXPEDU.ri.update, probs=0.975, na.rm=TRUE)

ri.lo <- c(CPI.ri.lo,GHSI.ri.lo,IGS.ri.lo,VAAG.ri.lo,GOVEXPEDU.ri.lo)
ri.med <- c(CPI.ri.med,GHSI.ri.med,IGS.ri.med,VAAG.ri.med,GOVEXPEDU.ri.med)
ri.up <- c(CPI.ri.up,GHSI.ri.up,IGS.ri.up,VAAG.ri.up,GOVEXPEDU.ri.up)

ri.out <- as.data.frame(cbind(ri.lo,ri.med,ri.up))
colnames(ri.out) <- c("ri.lo","ri.med","ri.up")
rownames(ri.out) <- attr(final.datPdam.imputed.orig.scale, "names")[predictors.sub]
ri.sort <- ri.out[order(ri.out[,2],decreasing=T),1:3]
ri.sort

pdam.final.phase.partial.pred.med <- pred.med
pdam.final.phase.partial.pred.up <- pred.up
pdam.final.phase.partial.pred.lo <- pred.lo
pdam.final.phase.partial.val.med <- val.med

# write outputs
dir.path <- paste(getwd(),"/out/",sep="")
for (o in 1:length(predictors.sub)) {
  var.name <- colnames(final.datPdam.imputed.orig.scale)[predictors.sub[o]]
  var.name.abbr <- toupper(sub("\\..*", "", var.name, ignore.case = T))
  var.title <- paste("pdam.final.phase.partial.",var.name.abbr,sep="")
  assign(var.title, 
         data.frame(val=pdam.final.phase.partial.val.med[,o],med=as.numeric(pdam.final.phase.partial.pred.med[,o]), 
                    up=as.numeric(pdam.final.phase.partial.pred.up[,o]), lo=as.numeric(pdam.final.phase.partial.pred.lo[,o])))
  get(var.title)
  write.table(get(var.title),paste(dir.path,paste(var.title,".csv",sep=""),sep=""), sep=",", row.names=F, col.names=T)
}

var.names <- rownames(ri.sort)
var.names.abbr <- toupper(sub("\\..*", "", var.names, ignore.case = T))
rel.infl.out <- data.frame(var=var.names.abbr, med=ri.sort$ri.med, up=ri.sort$ri.up, lo=ri.sort$ri.lo)
write.table(rel.infl.out,paste(dir.path,"pdam.final.phase.rel.infl.csv", sep=""), sep=",", row.names = F, col.names = T)




####################
####################
## response: r
####################
####################
head(final.datr.imputed.orig.scale)
hist(final.datr.imputed.orig.scale$r)

######################
# resampled BRT loop #
######################
biter <- 1000
eq.sp.points <- 100

## variable selection
colnames(final.datr.imputed.orig.scale)
pred.gdp.sub <- which(colnames(final.datr.imputed.orig.scale)=='gdp.sqr')
pred.cpi.sub <- which(colnames(final.datr.imputed.orig.scale)=='cpi.sqr')
pred.ghsi.sub <- which(colnames(final.datr.imputed.orig.scale)=='ghsi.sqr')
pred.igs.sub <- which(colnames(final.datr.imputed.orig.scale)=='igs.l')
pred.agrL.sub <- which(colnames(final.datr.imputed.orig.scale)=='agrL.l')
pred.VAag.sub <- which(colnames(final.datr.imputed.orig.scale)=='VAag.l')
pred.govexpedu.sub <- which(colnames(final.datr.imputed.orig.scale)=='govexpedu')
pred.stja.sub <- which(colnames(final.datr.imputed.orig.scale)=='stja.l')
predictors.sub <- c(pred.gdp.sub,pred.cpi.sub,pred.ghsi.sub,pred.igs.sub,pred.agrL.sub,
                    pred.VAag.sub,pred.govexpedu.sub,pred.stja.sub)
predictors.sub
head(final.datr.imputed.orig.scale[,predictors.sub])
resp.sub <- which(colnames(final.datr.imputed.orig.scale)=='r')
resp.sub
colnames(final.datr.imputed.orig.scale)[resp.sub]

# create storage arrays
val.arr <- pred.arr <- array(data = 0, dim = c(eq.sp.points, length(predictors.sub), biter),
                             dimnames=list(paste("x",1:eq.sp.points,sep=""),
                                           attr(final.datr.imputed.orig.scale, "names")[predictors.sub], paste("b",1:biter,sep="")))

# create storage vectors
D2.vec <- CV.cor.vec <- CV.cor.se.vec <- 
  GDP.ri <- CPI.ri <- GHSI.ri <- IGS.ri <- AGRL.ri <- VAAG.ri <- GOVEXPEDU.ri <- STJA.ri <- rep(0,biter)

for (b in 1:biter) {
  # resample data among countries
  resamp.sub <- sort(sample(x = 1:dim(final.datr.imputed.orig.scale)[1], size = dim(final.datr.imputed.orig.scale)[1], replace=TRUE))
  dat.resamp <- final.datr.imputed.orig.scale[resamp.sub,]
  
  # boosted regression tree
  brt.fit <- adaptive.gbm.step(dat.resamp, gbm.x = attr(dat.resamp, "names")[predictors.sub],
                               gbm.y = attr(dat.resamp, "names")[resp.sub], family="gaussian", max.trees=100000,
                               tolerance = 0.0001, learning.rate = 0.001, bag.fraction=0.75, tree.complexity = 2, silent=T, tolerance.method = "auto")
  summ.fit <- summary(brt.fit)
  
  # variable relative importance
  GDP.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[1]])])
  CPI.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[2]])])
  GHSI.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[3]])])
  IGS.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[4]])])
  AGRL.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[5]])])
  VAAG.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[6]])])
  GOVEXPEDU.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[7]])])
  STJA.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[8]])])
  
  D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / 
    brt.fit$cv.statistics$deviance.mean
  D2.vec[b] <- D2
  CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
  CV.cor.vec[b] <- CV.cor
  CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
  CV.cor.se.vec[b] <- CV.cor.se
  
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=length(predictors.sub))
  ## output average predictions
  for (p in 1:length(predictors.sub)) {
    RESP.val[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,2]
  }
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.fit$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.fit$var.names
  
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  print(b)
  
} # end b

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- 5
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] > (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  }
  print(k)
}

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)

val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median)

## plot
# order: GDP,CPI,GHSI,IGS,AGRL,VAAG,GOVEXPEDU,STJA
par(mfrow=c(2,4)) 

plot(val.med[,1],pred.med[,1],type="l",ylim=c(min(pred.lo[,1]),max(pred.up[,1])), lwd=2,
     ylab="(←lower) D:M ratio rate of change (higher→)", xlab="(←poorer) GDPpc (richer→)" )
lines(val.med[,1], pred.lo[,1], type="l", lty=2, col="red")
lines(val.med[,1], pred.up[,1], type="l", lty=2, col="red")

plot(val.med[,2],pred.med[,2],type="l",ylim=c(min(pred.lo[,2]),max(pred.up[,2])), lwd=2,
     ylab="(←lower) D:M ratio rate of change (higher→)", xlab="(←higher) corruption (lower→)" )
lines(val.med[,2], pred.lo[,2], type="l", lty=2, col="red")
lines(val.med[,2], pred.up[,2], type="l", lty=2, col="red")

plot(val.med[,3],pred.med[,3],type="l",ylim=c(min(pred.lo[,3]),max(pred.up[,3])), lwd=2,
     ylab="(←lower) D:M ratio rate of change (higher→)", xlab="(←lower) global health security (higher→)" )
lines(val.med[,3], pred.lo[,3], type="l", lty=3, col="red")
lines(val.med[,3], pred.up[,3], type="l", lty=3, col="red")

plot(val.med[,4],pred.med[,4],type="l",ylim=c(min(pred.lo[,4]),max(pred.up[,4])), lwd=2,
     ylab="(←lower) D:M ratio rate of change (higher→)", xlab="(←less) goods & services imports (more→)" )
lines(val.med[,4], pred.lo[,4], type="l", lty=3, col="red")
lines(val.med[,4], pred.up[,4], type="l", lty=3, col="red")

plot(val.med[,5],pred.med[,5],type="l",ylim=c(min(pred.lo[,5]),max(pred.up[,5])), lwd=2,
     ylab="(←lower) D:M ratio rate of change (higher→)", xlab="(←less) agricultural land (more→)" )
lines(val.med[,5], pred.lo[,5], type="l", lty=2, col="red")
lines(val.med[,5], pred.up[,5], type="l", lty=2, col="red")

plot(val.med[,6],pred.med[,6],type="l",ylim=c(min(pred.lo[,6]),max(pred.up[,6])), lwd=2,
     ylab="(←lower) D:M ratio rate of change (higher→)", xlab="(←less) value added primary production (more→)" )
lines(val.med[,6], pred.lo[,6], type="l", lty=2, col="red")
lines(val.med[,6], pred.up[,6], type="l", lty=2, col="red")

plot(val.med[,7],pred.med[,7],type="l",ylim=c(min(pred.lo[,7]),max(pred.up[,7])), lwd=2,
     ylab="(←lower) D:M ratio rate of change (higher→)", xlab="(←less) education investment (more→)" )
lines(val.med[,7], pred.lo[,7], type="l", lty=2, col="red")
lines(val.med[,7], pred.up[,7], type="l", lty=2, col="red")

plot(val.med[,8],pred.med[,8],type="l",ylim=c(min(pred.lo[,8]),max(pred.up[,8])), lwd=2,
     ylab="(←lower) D:M ratio rate of change (higher→)", xlab="(←fewer) scientific & tech journal articles (more→)" )
lines(val.med[,8], pred.lo[,8], type="l", lty=2, col="red")
lines(val.med[,8], pred.up[,8], type="l", lty=2, col="red")

par(mfrow=c(1,1)) 

# kappa method for output vectors
D2.update <- D2.vec[1:biter]
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# order: GDP,CPI,GHSI,IGS,AGRL,VAAG,GOVEXPEDU,STJA
GDP.ri.update <- GDP.ri[1:biter]
CPI.ri.update <- CPI.ri[1:biter]
GHSI.ri.update <- GHSI.ri[1:biter]
IGS.ri.update <- IGS.ri[1:biter]
AGRL.ri.update <- AGRL.ri[1:biter]
VAAG.ri.update <- VAAG.ri[1:biter]
GOVEXPEDU.ri.update <- GOVEXPEDU.ri[1:biter]
STJA.ri.update <- STJA.ri[1:biter]

for (k in 1:kappa.n) {
  D2.mean <- mean(D2.update, na.rm=T); D2.sd <- sd(D2.update, na.rm=T)
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  GDP.mean <- mean(GDP.ri.update, na.rm=T); GDP.sd <- sd(GDP.ri.update, na.rm=T)
  CPI.mean <- mean(CPI.ri.update, na.rm=T); CPI.sd <- sd(CPI.ri.update, na.rm=T)
  GHSI.mean <- mean(GHSI.ri.update, na.rm=T); GHSI.sd <- sd(GHSI.ri.update, na.rm=T)
  IGS.mean <- mean(IGS.ri.update, na.rm=T); IGS.sd <- sd(IGS.ri.update, na.rm=T)
  AGRL.mean <- mean(AGRL.ri.update, na.rm=T); AGRL.sd <- sd(AGRL.ri.update, na.rm=T)
  VAAG.mean <- mean(VAAG.ri.update, na.rm=T); VAAG.sd <- sd(VAAG.ri.update, na.rm=T)
  GOVEXPEDU.mean <- mean(GOVEXPEDU.ri.update, na.rm=T); GOVEXPEDU.sd <- sd(GOVEXPEDU.ri.update, na.rm=T)
  STJA.mean <- mean(STJA.ri.update, na.rm=T); STJA.sd <- sd(STJA.ri.update, na.rm=T)
  
  for (u in 1:biter) {
    D2.update[u] <- ifelse((D2.update[u] < (D2.mean-kappa*D2.sd) | D2.update[u] > (D2.mean+kappa*D2.sd)), NA, D2.update[u])
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] > (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] > (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    
    # order: GDP,CPI,GHSI,IGS,AGRL,VAAG,GOVEXPEDU,STJA
    GDP.ri.update[u] <- ifelse((GDP.ri.update[u] < (GDP.mean-kappa*GDP.sd) | GDP.ri.update[u] > (GDP.mean+kappa*GDP.sd)), NA, GDP.ri.update[u])
    CPI.ri.update[u] <- ifelse((CPI.ri.update[u] < (CPI.mean-kappa*CPI.sd) | CPI.ri.update[u] > (CPI.mean+kappa*CPI.sd)), NA, CPI.ri.update[u])
    GHSI.ri.update[u] <- ifelse((GHSI.ri.update[u] < (GHSI.mean-kappa*GHSI.sd) | GHSI.ri.update[u] > (GHSI.mean+kappa*GHSI.sd)), NA, GHSI.ri.update[u])
    IGS.ri.update[u] <- ifelse((IGS.ri.update[u] < (IGS.mean-kappa*IGS.sd) | IGS.ri.update[u] > (IGS.mean+kappa*IGS.sd)), NA, IGS.ri.update[u])
    AGRL.ri.update[u] <- ifelse((AGRL.ri.update[u] < (AGRL.mean-kappa*AGRL.sd) | AGRL.ri.update[u] > (AGRL.mean+kappa*AGRL.sd)), NA, AGRL.ri.update[u])
    VAAG.ri.update[u] <- ifelse((VAAG.ri.update[u] < (VAAG.mean-kappa*VAAG.sd) | VAAG.ri.update[u] > (VAAG.mean+kappa*VAAG.sd)), NA, VAAG.ri.update[u])
    GOVEXPEDU.ri.update[u] <- ifelse((GOVEXPEDU.ri.update[u] < (GOVEXPEDU.mean-kappa*GOVEXPEDU.sd) | GOVEXPEDU.ri.update[u] > (GOVEXPEDU.mean+kappa*GOVEXPEDU.sd)), NA, GOVEXPEDU.ri.update[u])
    STJA.ri.update[u] <- ifelse((STJA.ri.update[u] < (STJA.mean-kappa*STJA.sd) | STJA.ri.update[u] > (STJA.mean+kappa*STJA.sd)), NA, STJA.ri.update[u])
  }
  print(k)
}

D2.med <- median(D2.update, na.rm=TRUE)
D2.lo <- quantile(D2.update, probs=0.025, na.rm=TRUE)
D2.up <- quantile(D2.update, probs=0.975, na.rm=TRUE)
print(c(D2.lo,D2.med,D2.up))

CV.cor.med <- median(CV.cor.update, na.rm=TRUE)
CV.cor.lo <- quantile(CV.cor.update, probs=0.025, na.rm=TRUE)
CV.cor.up <- quantile(CV.cor.update, probs=0.975, na.rm=TRUE)
print(c(CV.cor.lo,CV.cor.med,CV.cor.up))

# order: GDP,CPI,GHSI,IGS,AGRL,VAAG,GOVEXPEDU,STJA
GDP.ri.lo <- quantile(GDP.ri.update, probs=0.025, na.rm=TRUE)
GDP.ri.med <- median(GDP.ri.update, na.rm=TRUE)
GDP.ri.up <- quantile(GDP.ri.update, probs=0.975, na.rm=TRUE)

CPI.ri.lo <- quantile(CPI.ri.update, probs=0.025, na.rm=TRUE)
CPI.ri.med <- median(CPI.ri.update, na.rm=TRUE)
CPI.ri.up <- quantile(CPI.ri.update, probs=0.975, na.rm=TRUE)

GHSI.ri.lo <- quantile(GHSI.ri.update, probs=0.025, na.rm=TRUE)
GHSI.ri.med <- median(GHSI.ri.update, na.rm=TRUE)
GHSI.ri.up <- quantile(GHSI.ri.update, probs=0.975, na.rm=TRUE)

IGS.ri.lo <- quantile(IGS.ri.update, probs=0.025, na.rm=TRUE)
IGS.ri.med <- median(IGS.ri.update, na.rm=TRUE)
IGS.ri.up <- quantile(IGS.ri.update, probs=0.975, na.rm=TRUE)

AGRL.ri.lo <- quantile(AGRL.ri.update, probs=0.025, na.rm=TRUE)
AGRL.ri.med <- median(AGRL.ri.update, na.rm=TRUE)
AGRL.ri.up <- quantile(AGRL.ri.update, probs=0.975, na.rm=TRUE)

VAAG.ri.lo <- quantile(VAAG.ri.update, probs=0.025, na.rm=TRUE)
VAAG.ri.med <- median(VAAG.ri.update, na.rm=TRUE)
VAAG.ri.up <- quantile(VAAG.ri.update, probs=0.975, na.rm=TRUE)

GOVEXPEDU.ri.lo <- quantile(GOVEXPEDU.ri.update, probs=0.025, na.rm=TRUE)
GOVEXPEDU.ri.med <- median(GOVEXPEDU.ri.update, na.rm=TRUE)
GOVEXPEDU.ri.up <- quantile(GOVEXPEDU.ri.update, probs=0.975, na.rm=TRUE)

STJA.ri.lo <- quantile(STJA.ri.update, probs=0.025, na.rm=TRUE)
STJA.ri.med <- median(STJA.ri.update, na.rm=TRUE)
STJA.ri.up <- quantile(STJA.ri.update, probs=0.975, na.rm=TRUE)

# order: GDP,CPI,GHSI,IGS,AGRL,VAAG,GOVEXPEDU,STJA
ri.lo <- c(GDP.ri.lo,CPI.ri.lo,GHSI.ri.lo,IGS.ri.lo,AGRL.ri.lo,VAAG.ri.lo,GOVEXPEDU.ri.lo,STJA.ri.lo)
ri.med <- c(GDP.ri.med,CPI.ri.med,GHSI.ri.med,IGS.ri.med,AGRL.ri.med,VAAG.ri.med,GOVEXPEDU.ri.med,STJA.ri.med)
ri.up <- c(GDP.ri.up,CPI.ri.up,GHSI.ri.up,IGS.ri.up,AGRL.ri.up,VAAG.ri.up,GOVEXPEDU.ri.up,STJA.ri.up)

ri.out <- as.data.frame(cbind(ri.lo,ri.med,ri.up))
colnames(ri.out) <- c("ri.lo","ri.med","ri.up")
rownames(ri.out) <- attr(final.datr.imputed.orig.scale, "names")[predictors.sub]
ri.sort <- ri.out[order(ri.out[,2],decreasing=T),1:3]
ri.sort

r.partial.pred.med <- pred.med
r.partial.pred.up <- pred.up
r.partial.pred.lo <- pred.lo
r.partial.val.med <- val.med

# write outputs
dir.path <- paste(getwd(),"/out/",sep="")
for (o in 1:length(predictors.sub)) {
  var.name <- colnames(final.datr.imputed.orig.scale)[predictors.sub[o]]
  var.name.abbr <- toupper(sub("\\..*", "", var.name, ignore.case = T))
  var.title <- paste("r.partial.",var.name.abbr,sep="")
  assign(var.title, 
         data.frame(val=r.partial.val.med[,o],med=as.numeric(r.partial.pred.med[,o]), 
                    up=as.numeric(r.partial.pred.up[,o]), lo=as.numeric(r.partial.pred.lo[,o])))
  get(var.title)
  write.table(get(var.title),paste(dir.path,paste(var.title,".csv",sep=""),sep=""), sep=",", row.names=F, col.names=T)
}

var.names <- rownames(ri.sort)
var.names.abbr <- toupper(sub("\\..*", "", var.names, ignore.case = T))
rel.infl.out <- data.frame(var=var.names.abbr, med=ri.sort$ri.med, up=ri.sort$ri.up, lo=ri.sort$ri.lo)
write.table(rel.infl.out,paste(dir.path,"r.rel.infl.csv", sep=""), sep=",", row.names = F, col.names = T)


## PHASE 1 - wealth/capacity
# gdp, cpi, ghsi, igs
colnames(final.datr.imputed.orig.scale)
pred.gdp.sub <- which(colnames(final.datr.imputed.orig.scale)=='gdp.sqr')
pred.cpi.sub <- which(colnames(final.datr.imputed.orig.scale)=='cpi.sqr')
pred.ghsi.sub <- which(colnames(final.datr.imputed.orig.scale)=='ghsi.sqr')
pred.igs.sub <- which(colnames(final.datr.imputed.orig.scale)=='igs.l')
predictors.sub <- c(pred.gdp.sub,pred.cpi.sub,pred.ghsi.sub,pred.igs.sub)
resp.sub <- which(colnames(final.datr.imputed.orig.scale)=='r')
head(final.datr.imputed.orig.scale[,predictors.sub])

brt.fit5 <- adaptive.gbm.step(final.datr.imputed.orig.scale, gbm.x = attr(final.datr.imputed.orig.scale, "names")[predictors.sub],
                              gbm.y = attr(final.datr.imputed.orig.scale, "names")[resp.sub],
                              family="gaussian", max.trees=100000, tolerance = 0.00003, learning.rate = 0.00001, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit5)
D2 <- 100 * (brt.fit5$cv.statistics$deviance.mean - brt.fit5$self.statistics$mean.resid) / brt.fit5$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit5)
gbm.plot.fits(brt.fit5)

brt5.CV.cor <- 100 * brt.fit5$cv.statistics$correlation.mean
brt5.CV.cor
brt5.CV.cor.se <- 100 * brt.fit5$cv.statistics$correlation.se
brt5.CV.cor.se
print(c(brt5.CV.cor, brt5.CV.cor.se))


# GLMM
# recode continent to greater region to increase per-category sample sizes
table(final.datr.imputed.orig.scale$cont)
final.datr.imputed.orig.scale$reg2 <- ifelse(final.datr.imputed.orig.scale$cont == "NAM" | final.datr.imputed.orig.scale$cont == "CAR", "NAMCAR", final.datr.imputed.orig.scale$cont)
final.datr.imputed.orig.scale$reg2 <- ifelse(final.datr.imputed.orig.scale$cont == "EUR" | final.datr.imputed.orig.scale$cont == "ME", "EURME", final.datr.imputed.orig.scale$reg2)
final.datr.imputed.orig.scale$reg2 <- ifelse(final.datr.imputed.orig.scale$cont == "ASIA" | final.datr.imputed.orig.scale$cont == "OC", "ASIAOC", final.datr.imputed.orig.scale$reg2)
table(final.datr.imputed.orig.scale$reg2)

# model set
colnames(final.datr.imputed.orig.scale)
vars <- c("gdp.sqr","cpi.sqr","ghsi.sqr","igs.l")    
vars4comb <- paste("r~", apply(combn(vars,4),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars3comb <- paste("r~", apply(combn(vars,3),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars2comb <- paste("r~", apply(combn(vars,2),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars1comb <- paste("r~", vars, "+(1|reg2)", sep="")
interconly <- "r~1+(1|reg2)"
mod.vec <- c(vars4comb,vars3comb,vars2comb,vars1comb,interconly)
mod.vec  

length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.datr.imputed.orig.scale, na.action=na.omit, REML=FALSE, verbose=F)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- AICc(fit)
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r2_nakagawa(fit)$R2_marginal # marginal R-squared
  Rc[i] <- 100*r2_nakagawa(fit)$R2_conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable1 <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable1) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable1) <- as.character(mod.vec)
summary.table1 <- sumtable1[order(sumtable1[,5],decreasing=F),]
summary.table1

fitsat <- lmer(as.formula(mod.vec[1]), data=final.datr.imputed.orig.scale, na.action=na.omit, REML=FALSE)
check_model(fitsat, detrend=T)
plot_model(fitsat, show.values=T, vline.color = "purple")
vif(fitsat)

## gdp collinear; remove
# model set
colnames(final.datr.imputed.orig.scale)
vars <- c("cpi.sqr","ghsi.sqr","igs.l")    
vars3comb <- paste("r~", apply(combn(vars,3),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars2comb <- paste("r~", apply(combn(vars,2),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars1comb <- paste("r~", vars, "+(1|reg2)", sep="")
interconly <- "r~1+(1|reg2)"
mod.vec <- c(vars3comb,vars2comb,vars1comb,interconly)
mod.vec  

length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.datr.imputed.orig.scale, na.action=na.omit, REML=FALSE, verbose=F)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- AICc(fit)
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r2_nakagawa(fit)$R2_marginal # marginal R-squared
  Rc[i] <- 100*r2_nakagawa(fit)$R2_conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable1 <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable1) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable1) <- as.character(mod.vec)
summary.table1 <- sumtable1[order(sumtable1[,5],decreasing=F),]
summary.table1

fitsat <- lmer(as.formula(mod.vec[1]), data=final.datr.imputed.orig.scale, na.action=na.omit, REML=FALSE)
check_model(fitsat, detrend=T)
plot_model(fitsat, show.values=T, vline.color = "purple")
vif(fitsat)


## PHASE 2 - primary production
# agrL, VAag
colnames(final.datr.imputed.orig.scale)
pred.agrl.sub <- which(colnames(final.datr.imputed.orig.scale)=='agrL.l')
pred.vaag.sub <- which(colnames(final.datr.imputed.orig.scale)=='VAag.l')
predictors.sub <- c(pred.agrl.sub,pred.vaag.sub)
resp.sub <- which(colnames(final.datr.imputed.orig.scale)=='r')

brt.fit1 <- adaptive.gbm.step(final.datr.imputed.orig.scale, gbm.x = attr(final.datr.imputed.orig.scale, "names")[predictors.sub],
                              gbm.y = attr(final.datr.imputed.orig.scale, "names")[resp.sub], 
                              family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.00001, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit1)
D2 <- 100 * (brt.fit1$cv.statistics$deviance.mean - brt.fit1$self.statistics$mean.resid) / brt.fit1$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit1)
gbm.plot.fits(brt.fit1)

brt1.CV.cor <- 100 * brt.fit1$cv.statistics$correlation.mean
brt1.CV.cor
brt1.CV.cor.se <- 100 * brt.fit1$cv.statistics$correlation.se
brt1.CV.cor.se
print(c(brt1.CV.cor, brt1.CV.cor.se))

# model set
vars <- c("agrL.l","VAag.l")    
vars2comb <- paste("r~", apply(combn(vars,2),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars1comb <- paste("r~", vars, "+(1|reg2)", sep="")
interconly <- "r~1+(1|reg2)"
mod.vec <- c(vars2comb,vars1comb,interconly)
mod.vec  

length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.datr.imputed.orig.scale, na.action=na.omit, REML=FALSE)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- AICc(fit)
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r2_nakagawa(fit)$R2_marginal # marginal R-squared
  Rc[i] <- 100*r2_nakagawa(fit)$R2_conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable2 <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable2) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable2) <- as.character(mod.vec)
summary.table2 <- sumtable2[order(sumtable2[,5],decreasing=F),]
summary.table2

fitsat <- lmer(as.formula(mod.vec[1]), data=final.datr.imputed.orig.scale, na.action=na.omit, REML=FALSE)
check_model(fitsat, detrend=T)
plot_model(fitsat, show.values=T, vline.color = "purple")


## PHASE 3 - education/research
# govexpedu, stja
colnames(final.datr.imputed.orig.scale)
pred.govexpedu.sub <- which(colnames(final.datr.imputed.orig.scale)=='govexpedu')
pred.stja.l.sub <- which(colnames(final.datr.imputed.orig.scale)=='stja.l')
predictors.sub <- c(pred.govexpedu.sub,pred.stja.l.sub)
resp.sub <- which(colnames(final.datr.imputed.orig.scale)=='r')

brt.fit1 <- adaptive.gbm.step(final.datr.imputed.orig.scale, gbm.x = attr(final.datr.imputed.orig.scale, "names")[predictors.sub],
                              gbm.y = attr(final.datr.imputed.orig.scale, "names")[resp.sub],
                              family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.00001, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit1)
D2 <- 100 * (brt.fit1$cv.statistics$deviance.mean - brt.fit1$self.statistics$mean.resid) / brt.fit1$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit1)
gbm.plot.fits(brt.fit1)

brt1.CV.cor <- 100 * brt.fit1$cv.statistics$correlation.mean
brt1.CV.cor
brt1.CV.cor.se <- 100 * brt.fit1$cv.statistics$correlation.se
brt1.CV.cor.se
print(c(brt1.CV.cor, brt1.CV.cor.se))


# GLMM
# model set
colnames(final.datr.imputed.orig.scale)
vars <- c("govexpedu","stja.l")    
vars2comb <- paste("r~", apply(combn(vars,2),2,paste,collapse='+'),"+(1|reg2)", sep="")
vars1comb <- paste("r~", vars, "+(1|reg2)", sep="")
interconly <- "r~1+(1|reg2)"
mod.vec <- c(vars2comb,vars1comb,interconly)
mod.vec  

length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.datr.imputed.orig.scale, na.action=na.omit, REML=FALSE)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- AIC(fit)
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r2_nakagawa(fit)$R2_marginal # marginal R-squared
  Rc[i] <- 100*r2_nakagawa(fit)$R2_conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable3 <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable3) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable3) <- as.character(mod.vec)
summary.table3 <- sumtable3[order(sumtable3[,5],decreasing=F),]
summary.table3

fitsat <- lmer(as.formula(mod.vec[1]), data=final.datr.imputed.orig.scale, na.action=na.omit, REML=FALSE)
check_model(fitsat, detrend=T)
plot_model(fitsat, show.values=T, vline.color = "purple")



## COMBINE PHASES
# cpi, ghsi, igs, agrL, govexpedu
# model set
colnames(final.datr.imputed.orig.scale)
vars <- c("cpi.sqr", "ghsi.sqr", "igs.l", "agrL.l", "govexpedu")

vars5comb <- paste0(
  "r~",
  paste(vars, collapse = "+"),
  "+(1|reg2)"
)

vars4comb <- paste(
  "r~",
  apply(combn(vars, 4), 2, paste, collapse = "+"),
  "+(1|reg2)",
  sep = ""
)

vars3comb <- paste(
  "r~",
  apply(combn(vars, 3), 2, paste, collapse = "+"),
  "+(1|reg2)",
  sep = ""
)

vars2comb <- paste(
  "r~",
  apply(combn(vars, 2), 2, paste, collapse = "+"),
  "+(1|reg2)",
  sep = ""
)

vars1comb <- paste0(
  "r~",
  vars,
  "+(1|reg2)"
)

interconly <- "r~1+(1|reg2)"

mod.vec <- unique(c(
  vars5comb,
  vars4comb,
  vars3comb,
  vars2comb,
  vars1comb,
  interconly
))

stopifnot(
  length(mod.vec) == 2^length(vars),
  length(mod.vec) == length(unique(mod.vec))
)

length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=final.datr.imputed.orig.scale, na.action=na.omit, REML=FALSE)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- AICc(fit)
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r2_nakagawa(fit)$R2_marginal # marginal R-squared
  Rc[i] <- 100*r2_nakagawa(fit)$R2_conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtableF <- data.frame(mod.num,k.vec,round(LL.vec,3),round(AICc.vec,3),round(dAICc,3),round(wAICc,3),
                        round(BIC.vec,3),round(dBIC,3),round(wBIC,4),round(Rm,1),round(Rc,1))
colnames(sumtableF) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtableF) <- as.character(mod.vec)
summary.tableF <- sumtableF[order(sumtableF[,5],decreasing=F),]
summary.tableF

fitsat <- lmer(as.formula(mod.vec[1]), data=final.datr.imputed.orig.scale, na.action=na.omit, REML=FALSE)
vif(fitsat)
check_model(fitsat, detrend=F)
plot_model(fitsat, type="est", sort.est=T)


## BRT
# cpi, ghsi, igs, agrL, govexpedu
colnames(final.datr.imputed.orig.scale)
pred.cpi.sub <- which(colnames((final.datr.imputed.orig.scale))=='cpi.sqr')
pred.ghsi.sub <- which(colnames((final.datr.imputed.orig.scale))=='ghsi.sqr')
pred.igs.sub <- which(colnames(final.datr.imputed.orig.scale)=='igs.l')
pred.agrL.sub <- which(colnames(final.datr.imputed.orig.scale)=='agrL.l')
pred.govexpedu.sub <- which(colnames(final.datr.imputed.orig.scale)=='govexpedu')
predictors.sub <- c(pred.cpi.sub,pred.ghsi.sub,pred.igs.sub,pred.agrL.sub,pred.govexpedu.sub)
resp.sub <- which(colnames(final.datr.imputed.orig.scale)=='r')

brt.fitF <- adaptive.gbm.step(final.datr.imputed.orig.scale, gbm.x = attr(final.datr.imputed.orig.scale, "names")[predictors.sub],
                              gbm.y = attr(final.datr.imputed.orig.scale, "names")[resp.sub], family="gaussian", 
                              max.trees=100000, tolerance = 0.00001, learning.rate = 0.00003, 
                              bag.fraction=0.75, tree.complexity = 2)
summary(brt.fitF)
D2 <- 100 * (brt.fitF$cv.statistics$deviance.mean - brt.fitF$self.statistics$mean.resid) / brt.fitF$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fitF)
gbm.plot.fits(brt.fitF)

brtF.CV.cor <- 100 * brt.fitF$cv.statistics$correlation.mean
brtF.CV.cor
brtF.CV.cor.se <- 100 * brt.fitF$cv.statistics$correlation.se
brtF.CV.cor.se
print(c(brtF.CV.cor, brtF.CV.cor.se))


# resampled BRT loop
# cpi, ghsi, igs, agrL, govexpedu
biter <- 200
eq.sp.points <- 100

colnames(final.datr.imputed.orig.scale)
pred.cpi.sub <- which(colnames((final.datr.imputed.orig.scale))=='cpi.sqr')
pred.ghsi.sub <- which(colnames((final.datr.imputed.orig.scale))=='ghsi.sqr')
pred.igs.sub <- which(colnames(final.datr.imputed.orig.scale)=='igs.l')
pred.agrL.sub <- which(colnames(final.datr.imputed.orig.scale)=='agrL.l')
pred.govexpedu.sub <- which(colnames(final.datr.imputed.orig.scale)=='govexpedu')
predictors.sub <- c(pred.cpi.sub,pred.ghsi.sub,pred.igs.sub,pred.agrL.sub,pred.govexpedu.sub)
colnames(final.datr.imputed.orig.scale[,predictors.sub])
resp.sub <- which(colnames(final.datr.imputed.orig.scale)=='r')

# create storage arrays
val.arr <- pred.arr <- array(data = NA, dim = c(eq.sp.points, length(predictors.sub), biter),
                             dimnames=list(paste("x",1:eq.sp.points,sep=""), 
                                           attr(final.datr.imputed.orig.scale, "names")[predictors.sub], 
                                           paste("b",1:biter,sep="")))

# create storage vectors
D2.vec <- CV.cor.vec <- CV.cor.se.vec <- CPI.ri <- GHSI.ri <- IGS.ri <- 
  AGRL.ri <- GOVEXPEDU.ri <- rep(NA,biter)

for (b in 1:biter) {
  # resample data among countries
  resamp.sub <- sort(sample(x = 1:dim(final.datr.imputed.orig.scale)[1], 
                            size = dim(final.datr.imputed.orig.scale)[1], replace=TRUE))
  dat.resamp <- final.datr.imputed.orig.scale[resamp.sub,]
  
  # boosted regression tree
  brt.fit <- adaptive.gbm.step(dat.resamp, gbm.x = attr(dat.resamp, "names")[predictors.sub],
                               gbm.y = attr(dat.resamp, "names")[resp.sub], family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.0001, bag.fraction=0.75, tree.complexity = 2, silent=T, tolerance.method = "auto")
  summ.fit <- summary(brt.fit)
  length(summ.fit[[1]])
  
  # variable relative importance
  # cpi, ghsi, igs, agrL, govexpedu
  CPI.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[1]])])
  GHSI.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[2]])])
  IGS.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[3]])])
  AGRL.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[4]])])
  GOVEXPEDU.ri[b] <- as.numeric(summ.fit$rel.inf[which(rownames(summ.fit) == colnames(dat.resamp)[predictors.sub[5]])])

  D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
  D2.vec[b] <- D2
  CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
  CV.cor.vec[b] <- CV.cor
  CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
  CV.cor.se.vec[b] <- CV.cor.se
  
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=length(predictors.sub))
  ## output average predictions
  for (p in 1:length(predictors.sub)) {
    RESP.val[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,2]
  }
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.fit$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.fit$var.names
  
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  print(b)
  
} # end b

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- 5
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  }
  print(k)
}

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)

val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

par(mfrow=c(2,3)) 
# cpi, ghsi, igs, agrL, govexpedu
plot(val.med[,1],pred.med[,1],type="l",ylim=c(min(pred.lo[,1]),max(pred.up[,1])),
     lwd=2, ylab="(←lower) dam:mgm r (higher→)",  xlab="(←higher) corruption (lower→)" )
lines(val.med[,1], pred.lo[,1], type="l", lty=2, col="red")
lines(val.med[,1], pred.up[,1], type="l", lty=2, col="red")

plot(val.med[,2],pred.med[,2],type="l",ylim=c(min(pred.lo[,2]),max(pred.up[,2])),
     lwd=2, ylab="(←lower) dam:mgm r (higher→)",  xlab="(←lower) global health security index (higher→)" )
lines(val.med[,2], pred.lo[,2], type="l", lty=2, col="red")
lines(val.med[,2], pred.up[,2], type="l", lty=2, col="red")

plot(val.med[,3],pred.med[,3],type="l",ylim=c(min(pred.lo[,3]),max(pred.up[,3])),
     lwd=2, ylab="(←lower) dam:mgm r (higher→)",  xlab="(←lower) per capita imports goods & services (higher→)" )
lines(val.med[,3], pred.lo[,3], type="l", lty=2, col="red")
lines(val.med[,3], pred.up[,3], type="l", lty=2, col="red")

plot(val.med[,4],pred.med[,4],type="l",ylim=c(min(pred.lo[,4]),max(pred.up[,4])),
     lwd=3, ylab="(←lower) dam:mgm r (higher→)", xlab="(←less) % agricultural land (more→)")
lines(val.med[,4], pred.lo[,4], type="l", lty=2, col="red")
lines(val.med[,4], pred.up[,4], type="l", lty=2, col="red")

plot(val.med[,5],pred.med[,5],type="l",ylim=c(min(pred.lo[,5]),max(pred.up[,5])),
     lwd=3, ylab="(←lower) dam:mgm r (higher→)", xlab="(←fewer) gov education expenditure (more→)")
lines(val.med[,5], pred.lo[,5], type="l", lty=2, col="red")
lines(val.med[,5], pred.up[,5], type="l", lty=2, col="red")

par(mfrow=c(1,1)) 

# kappa method for output vectors
D2.update <- D2.vec[1:biter]
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# igs, agrL, govexpedu
CPI.ri.update <- CPI.ri[1:biter]
GHSI.ri.update <- GHSI.ri[1:biter]
IGS.ri.update <- IGS.ri[1:biter]
AGRL.ri.update <- AGRL.ri[1:biter]
GOVEXPEDU.ri.update <- GOVEXPEDU.ri[1:biter]

for (k in 1:kappa.n) {
  D2.mean <- mean(D2.update, na.rm=T); D2.sd <- sd(D2.update, na.rm=T)
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  CPI.mean <- mean(CPI.ri.update, na.rm=T); CPI.sd <- sd(CPI.ri.update, na.rm=T)
  GHSI.mean <- mean(GHSI.ri.update, na.rm=T); GHSI.sd <- sd(GHSI.ri.update, na.rm=T)
  IGS.mean <- mean(IGS.ri.update, na.rm=T); IGS.sd <- sd(IGS.ri.update, na.rm=T)
  AGRL.mean <- mean(AGRL.ri.update, na.rm=T); AGRL.sd <- sd(AGRL.ri.update, na.rm=T)
  GOVEXPEDU.mean <- mean(GOVEXPEDU.ri.update, na.rm=T); GOVEXPEDU.sd <- sd(GOVEXPEDU.ri.update, na.rm=T)
  
  for (u in 1:biter) {
    D2.update[u] <- ifelse((D2.update[u] < (D2.mean-kappa*D2.sd) | D2.update[u] > (D2.mean+kappa*D2.sd)), NA, D2.update[u])
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] > (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] > (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    
    CPI.ri.update[u] <- ifelse((CPI.ri.update[u] < (CPI.mean-kappa*CPI.sd) | CPI.ri.update[u] > (CPI.mean+kappa*CPI.sd)), NA, CPI.ri.update[u])
    GHSI.ri.update[u] <- ifelse((GHSI.ri.update[u] < (GHSI.mean-kappa*GHSI.sd) | GHSI.ri.update[u] > (GHSI.mean+kappa*GHSI.sd)), NA, GHSI.ri.update[u])
    IGS.ri.update[u] <- ifelse((IGS.ri.update[u] < (IGS.mean-kappa*IGS.sd) | IGS.ri.update[u] > (IGS.mean+kappa*IGS.sd)), NA, IGS.ri.update[u])
    AGRL.ri.update[u] <- ifelse((AGRL.ri.update[u] < (AGRL.mean-kappa*AGRL.sd) | AGRL.ri.update[u] > (AGRL.mean+kappa*AGRL.sd)), NA, AGRL.ri.update[u])
    GOVEXPEDU.ri.update[u] <- ifelse((GOVEXPEDU.ri.update[u] < (GOVEXPEDU.mean-kappa*GOVEXPEDU.sd) | GOVEXPEDU.ri.update[u] > (GOVEXPEDU.mean+kappa*GOVEXPEDU.sd)), NA, GOVEXPEDU.ri.update[u])
  }
  
  print(k)
}

D2.med <- median(D2.update, na.rm=TRUE)
D2.lo <- quantile(D2.update, probs=0.025, na.rm=TRUE)
D2.up <- quantile(D2.update, probs=0.975, na.rm=TRUE)
print(c(D2.lo,D2.med,D2.up))

CV.cor.med <- median(CV.cor.update, na.rm=TRUE)
CV.cor.lo <- quantile(CV.cor.update, probs=0.025, na.rm=TRUE)
CV.cor.up <- quantile(CV.cor.update, probs=0.975, na.rm=TRUE)
print(c(CV.cor.lo,CV.cor.med,CV.cor.up))

CPI.ri.lo <- quantile(CPI.ri.update, probs=0.025, na.rm=TRUE)
CPI.ri.med <- median(CPI.ri.update, na.rm=TRUE)
CPI.ri.up <- quantile(CPI.ri.update, probs=0.975, na.rm=TRUE)

GHSI.ri.lo <- quantile(GHSI.ri.update, probs=0.025, na.rm=TRUE)
GHSI.ri.med <- median(GHSI.ri.update, na.rm=TRUE)
GHSI.ri.up <- quantile(GHSI.ri.update, probs=0.975, na.rm=TRUE)

IGS.ri.lo <- quantile(IGS.ri.update, probs=0.025, na.rm=TRUE)
IGS.ri.med <- median(IGS.ri.update, na.rm=TRUE)
IGS.ri.up <- quantile(IGS.ri.update, probs=0.975, na.rm=TRUE)

AGRL.ri.lo <- quantile(AGRL.ri.update, probs=0.025, na.rm=TRUE)
AGRL.ri.med <- median(AGRL.ri.update, na.rm=TRUE)
AGRL.ri.up <- quantile(AGRL.ri.update, probs=0.975, na.rm=TRUE)

GOVEXPEDU.ri.lo <- quantile(GOVEXPEDU.ri.update, probs=0.025, na.rm=TRUE)
GOVEXPEDU.ri.med <- median(GOVEXPEDU.ri.update, na.rm=TRUE)
GOVEXPEDU.ri.up <- quantile(GOVEXPEDU.ri.update, probs=0.975, na.rm=TRUE)

ri.lo <- c(CPI.ri.lo,GHSI.ri.lo,IGS.ri.lo,AGRL.ri.lo,GOVEXPEDU.ri.lo)
ri.med <- c(CPI.ri.med,GHSI.ri.med,IGS.ri.med,AGRL.ri.med,GOVEXPEDU.ri.med)
ri.up <- c(CPI.ri.up,GHSI.ri.up,IGS.ri.up,AGRL.ri.up,GOVEXPEDU.ri.up)

ri.out <- as.data.frame(cbind(ri.lo,ri.med,ri.up))
colnames(ri.out) <- c("ri.lo","ri.med","ri.up")
rownames(ri.out) <- attr(final.datr.imputed.orig.scale, "names")[predictors.sub]
ri.sort <- ri.out[order(ri.out[,2],decreasing=T),1:3]
ri.sort

r.final.phase.partial.pred.med <- pred.med
r.final.phase.partial.pred.up <- pred.up
r.final.phase.partial.pred.lo <- pred.lo
r.final.phase.partial.val.med <- val.med

# write outputs
dir.path <- paste(getwd(),"/out/",sep="")
for (o in 1:length(predictors.sub)) {
  var.name <- colnames(final.datr.imputed.orig.scale)[predictors.sub[o]]
  var.name.abbr <- toupper(sub("\\..*", "", var.name, ignore.case = T))
  var.title <- paste("r.final.phase.partial.",var.name.abbr,sep="")
  assign(var.title, 
         data.frame(val=r.final.phase.partial.val.med[,o],med=as.numeric(r.final.phase.partial.pred.med[,o]), 
                    up=as.numeric(r.final.phase.partial.pred.up[,o]), lo=as.numeric(r.final.phase.partial.pred.lo[,o])))
  get(var.title)
  write.table(get(var.title),paste(dir.path,paste(var.title,".csv",sep=""),sep=""), sep=",", row.names=F, col.names=T)
}

var.names <- rownames(ri.sort)
var.names.abbr <- toupper(sub("\\..*", "", var.names, ignore.case = T))
rel.infl.out <- data.frame(var=var.names.abbr, med=ri.sort$ri.med, up=ri.sort$ri.up, lo=ri.sort$ri.lo)
write.table(rel.infl.out,paste(dir.path,"r.final.phase.rel.infl.csv", sep=""), sep=",", row.names = F, col.names = T)



#######################################################################################
## government expenditure in education as a function of economic indicators (IGS & GDP)
#######################################################################################
plot(final.dat.imputed.orig.scale$igs.l, final.dat.imputed.orig.scale$govexpedu, pch=19, 
     ylab="government education expenditure (%GDP)", xlab="imports of goods & services")
abline(lm(final.dat.imputed.orig.scale$govexpedu ~ final.dat.imputed.orig.scale$igs.l), lty=2, col="red")

plot(final.dat.imputed.orig.scale$gdp.sqr, final.dat.imputed.orig.scale$govexpedu, pch=19, 
     ylab="government education expenditure (%GDP)", xlab="per-capita GDP")
abline(lm(final.dat.imputed.orig.scale$govexpedu ~ final.dat.imputed.orig.scale$gdp.sqr), lty=2, col="red")
