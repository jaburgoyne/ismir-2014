### ANOVA analysis of MIREX 2013 results.
### Author: J. A. Burgoyne
### Date: November 2013

library(boot)     # For inv.logit().
library(car)      # For contr.Sum(), which is more readable than contr.sum().
library(multcomp) # For glht() and mcp().

###
### I: Set constants.
###

ALGOS <- c("CB3", "CB4", "CF2", "KO1", "KO2", "NG1",
           "NG2", "NMSD1", "NMSD2", "PP3", "PP4", "SB8")
VOCABS <- c("Root", "MajMin", "MajMinBass", "Sevenths", "SeventhsBass")
ANALYSES <- c(VOCABS, "Segmentation")
## Each analysis requires its own duration of evaluated
## symbols. Segmentation should use the duration of the entire song,
## which should be equivalent to DurRoot.
DURATIONS <- c(sapply(VOCABS, function (x) paste0("Dur", x)),
               list(Segmentation = "DurRoot"))
DATA.SETS <- c("MirexChord2009", "BillboardTest2012", "BillboardTest2013")

###
### II: Read Johan Pauwels's result tables.
###

## Reads the chord symbol recall for all algorithms on a
## particular vocabulary (or segmentation results for all
## algorithms). Must be run from within the results directory of a
## specific data set.
read.results <- function (analysis) {
  isSeg <- analysis == "Segmentation"
  cnames <- c()
  cclasses <- c()
  if (isSeg) {
    setwd("resultsSegmentation")
    cnames <- c("Song", "", "UnderSeg", "OverSeg")
    cclasses <- c("factor", "NULL", "numeric", "numeric")
  } else {
    setwd(paste0("resultsMirex", analysis))
    cnames <- c("Song", analysis, paste0("Dur", analysis), rep("", 6))
    cclasses <- c("factor", "numeric", "numeric", rep("NULL", 6))
  }
  results <- data.frame()
  for (algo in ALGOS) {
    algo.results <- read.table(paste0(algo, ".csv"),
                               skip = 2,
                               header = FALSE,
                               sep = ",",
                               dec = ".",
                               col.names = cnames, 
                               colClasses = cclasses)
    if (isSeg) {
      ## Replace over- and under-segmentation with the harmonic mean
      ## of their arithmetic inverses.
      algo.results$Segmentation <-
        1 / (0.5 * (1/(1-algo.results$UnderSeg) + 1/(1-algo.results$OverSeg)))
      algo.results$OverSeg <- c()
      algo.results$UnderSeg <- c()
    } else {
      ## Rescale RCO results back to [0, 1].
      algo.results[,2] <- algo.results[,2]/100
    }
    results <- rbind(results,
                     cbind(Algo = rep(algo, dim(algo.results)[1]),
                           algo.results))
  }
  setwd("..")
  return(results)
}

## Reads the CSR and segmentation data for all vocabularies.
read.all.results <- function (data.set) {
  setwd(data.set)
  results <- Reduce(function (x, y) merge(x, y, all=TRUE),
                    lapply(ANALYSES, read.results))
  ## Set the contrasts such that they represent differences from the
  ## grand average.
  contrasts(results$Song) <- contr.Sum(levels(results$Song))
  contrasts(results$Algo) <- contr.Sum(levels(results$Algo))
  setwd("..")
  return(results)
}

## List of results for each data set.
results <- sapply(DATA.SETS, read.all.results, simplify=FALSE)

###
### III: Fit logit models.
###

## Fits a logit model including song and algorithm. Requires a data
## frame with four columns: `p` for the probability being fit, `song`
## for the song, `algo` for the algorithm, and `dur` for the duration
## of the evaluated symbols used to compute `p`. Returns a glht object
## including pairwise comparisons for all algorithms.
comparisons <- function (dat) {
  ## ANOVA. We want logit models (the default) rather than probit
  ## models because they are more robust to over-dispersion (heavy
  ## tails).
  prefit <- glm(p ~ song + algo,
                data = dat,
                family = quasibinomial(link = "logit"),
                weights = dur)
  ## Reorder the levels of the algorithm factor in descending order of
  ## performance in order to facilitate pairwise comparisions. With
  ## sum-to-one contrasts, we can recover the value of the aliased
  ## last level from the arithmetic inverse of the sum of the
  ## coefficients for all other levels.
  M <- nlevels(dat$song)
  N <- nlevels(dat$algo)
  algo.coefs <- coef(prefit)[c(M + 1:(N-1), 1)]
  algo.coefs[N] <- -sum(algo.coefs[1:(N-1)])
  new.algo <-
    factor(dat$algo,
           levels = levels(dat$algo)[order(algo.coefs, decreasing=TRUE)])
  contrasts(new.algo) <- contr.Sum(levels(new.algo))
  fit <- glm(p ~ song + new.algo,
             data = dat,
             family = quasibinomial(link = "logit"),
             weights = dur)
  ## Compute the differences between all pairwise combinations of
  ## algorithms.
  tuk <- glht(fit, mcp(new.algo = "Tukey"))
}

## List of lists of pairwise comparisons for each data set and
## vocabulary or segmentation.
all.comparisons <-
  sapply(results,
         function (dat) {
           sapply(ANALYSES,
                  function (analysis) {
                    comparisons(data.frame(p = dat[,analysis],
                                           song = dat$Song,
                                           algo = dat$Algo,
                                           dur = dat[,DURATIONS[[analysis]]]))
                  },
                  simplify = FALSE)
         },
         simplify = FALSE)

###
### IV: Plotting functions
###

## Draft function for plotting results. Based on plot.glht, which is
## in turn based on plot.TukeyHSD.
plot.odds <- function(y, xlim, xlab, ylim, ...) {
    x <- confint(y)
    xi <- exp(x$confint)
    dimnames(xi)[[1]] <- sub("-", "/", dimnames(xi)[[1]])
    ### make sure one-sided intervals are drawn correctly
    xrange <- c(min(xi[,"lwr"]), max(xi[, "upr"]))
    if (!is.finite(xrange[1])) xrange[1] <- min(xi[,"Estimate"])
    if (!is.finite(xrange[2])) xrange[2] <- max(xi[,"Estimate"])
    yvals <- nrow(xi):1
    if (missing(xlim))
        xlim <- xrange
    if (missing(ylim))
        ylim <- c(0.5, nrow(xi) + 0.5)
    plot(c(xi[, "lwr"], xi[, "upr"]), rep.int(yvals, 2), 
         type = "n", axes = FALSE, xlab = "", ylab = "", 
         xlim = xlim, ylim = ylim, ...)
    axis(1, ...)
    axis(2, at = nrow(xi):1, labels = dimnames(xi)[[1]], 
         las = 1, ...)
    abline(h = yvals, lty = 1, lwd = 1, col = "lightgray")
    abline(v = 1, lty = 2, lwd = 1, ...)
    left <- xi[, "lwr"]
    left[!is.finite(left)] <- min(c(0, xlim[1] * 2))
    right <- xi[, "upr"]
    right[!is.finite(right)] <- max(c(0, xlim[2] * 2))
    segments(left, yvals, right, yvals, ...)
    points(xi[, "lwr"], yvals, pch = "(", ...)
    points(xi[, "upr"], yvals, pch = ")", ...)
    points(xi[, "Estimate"], yvals, pch = 20, ...)
    main <- list(...)$main
    if (is.null(main)) {
        if (attr(x, "type") == "adjusted") {
            main <- paste(format(100 * attr(x$confint, "conf.level"), 2), 
                          "% Family-Wise Confidence Level\n", sep = "")
        } else {
            main <- paste(format(100 * attr(x$confint, "conf.level"), 2),
                          "% Confidence Level\n", sep = "")
        }
    } else {
        main <- NULL ### main was already plotted in plot() via ...
    }
    if (missing(xlab))
          xlab <- "Odds ratio"
    title(main = main, xlab = xlab)
    box()
}

## Logistic probability plot
qqlogis <- function (y) {
  qqplot(qlogis(ppoints(y)), y,
         xlab="Theoretical Quantiles", ylab="Sample Quantiles")
  qqline(y, distribution=qlogis)
}
  
