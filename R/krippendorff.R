
krippendorff.within = function(data, dist)
{
    m = rowSums(! is.na(data))
    D.o = 0
    n.u = nrow(data)
    n.c = ncol(data)
    for (i in 1:n.u)
    {
        s = 0
        if (m[i] > 1)
        {
            for (j in 1:(n.c - 1))
                for (k in (j + 1):n.c)
                    s = s + dist(data[i, j], data[i, k])
            D.o = D.o + s / (m[i] - 1)
        }
    }
    D.o = D.o / sum(m[m > 1])
}

krippendorff.total = function(data, dist)
{
    m = rowSums(! is.na(data))
    D.e = 0
    n.u = nrow(data)
    n.c = ncol(data)
    for (i in 1:n.u)
        for (j in 1:n.c)
            for (k in 1:n.u)
                for (l in 1:n.c)
                    D.e = D.e + dist(data[i, j], data[k, l])
    n = sum(m)
    D.e = D.e / (2 * n * (n - 1))
}

krippendorff.control = function(confint, verbose, control)
{
	if (length(control) > 0)
	{
        nms = match.arg(names(control), c("bootit", "parallel", "type", "nodes"), several.ok = TRUE)
        control = control[nms]
    }
    bootit = control$bootit
    if (confint)
    {
        bootit = control$bootit
        if (is.null(bootit) || ! is.numeric(bootit) || length(bootit) > 1 || bootit != as.integer(bootit) || bootit < 100)
        {
            if (verbose)
                message("\nControl parameter 'bootit' must be a positive integer >= 100. Setting it to the default value of 1,000.")
            control$bootit = 1000
        }
        if (is.null(control$parallel) || ! is.logical(control$parallel) || length(control$parallel) > 1)
        {
            if (verbose)
                message("\nControl parameter 'parallel' must be a logical value. Setting it to the default value of TRUE.")
            control$parallel = TRUE
        }
        if (control$parallel)
        {
            if (requireNamespace("parallel", quietly = TRUE))
            {
                if (is.null(control$type) || length(control$type) > 1 || ! control$type %in% c("SOCK", "PVM", "MPI", "NWS"))
                {
                    if (verbose)
                        message("\nControl parameter 'type' must be \"SOCK\", \"PVM\", \"MPI\", or \"NWS\". Setting it to \"SOCK\".")
                    control$type = "SOCK"
                }
                nodes = control$nodes
                if (is.null(control$nodes) || ! is.numeric(nodes) || length(nodes) > 1 || nodes != as.integer(nodes) || nodes < 2)
                    stop("Control parameter 'nodes' must be a whole number greater than 1.")
            }
            else
            {
                if (verbose)
                    message("\nParallel computation requires package parallel. Setting control parameter 'parallel' to FALSE.")
                control$parallel = FALSE
                control$type = control$nodes = NULL 
            }
        }
    }
    else
        control$bootit = control$parallel = control$type = control$nodes = NULL
    control
}

krippendorff.bootstrap = function(object)
{
    boot.sample = numeric(object$control$bootit)
    data = vector("list", object$control$bootit)
    for (j in 1:object$control$bootit)
        data[[j]] = object$data[sample(nrow(object$data), replace = TRUE), ]
    if (! object$control$parallel)
    {
        if (object$verbose && requireNamespace("pbapply", quietly = TRUE))
        {
            gathered = pbapply::pblapply(data, krippendorff.within, dist = object$dist)
            for (j in 1:object$control$bootit)
            {
                D.o = gathered[[j]]
                boot.sample[j] = 1 - D.o / object$D.e
            }
        }
        else
        {
            for (j in 1:object$control$bootit)
            {
                D.o = krippendorff.within(data[[j]], dist = object$dist)
                boot.sample[j] = 1 - D.o / object$D.e
            }
        }
    }
    else
    {
        cl = parallel::makeCluster(object$control$nodes, object$control$type)
        if (object$verbose && requireNamespace("pbapply", quietly = TRUE))
            gathered = pbapply::pblapply(data, krippendorff.within, dist = object$dist, cl = cl)
        else
            gathered = parallel::clusterApplyLB(cl, data, krippendorff.within, dist = object$dist)
        parallel::stopCluster(cl)
        for (j in 1:object$control$bootit)
        {
            D.o = gathered[[j]]
            boot.sample[j] = 1 - D.o / object$D.e
        }
    }
    boot.sample
}

#' Apply Krippendorff's Alpha.
#'
#' @details This is the package's flagship function. It applies the Krippendorff's Alpha methodology for nominal, ordinal, interval, or ratio levels of measurement, and, if desired, produces confidence intervals. Parallel computing is supported, when applicable.
#'
#' If the level of measurement is nominal, the discrete metric (\code{\link{nominal.dist}}) is employed by default. If the level of measurement is interval or ordinal, the squared-difference distance function (\code{\link{interval.dist}}) is employed by default. (For the ordinal level of measurement, using the squared-difference distance function may be inappropriate, in which case the user should supply his/her own distance function.) If the level of measurement is ratio, a ratio distance function (\code{\link{ratio.dist}}) is applied. Alternatively, the user may supply his/her own distance function. Said function must handle \code{NA}'s gracefully; see the above mentioned built-in distance functions for examples.
#'
#' If argument \code{confint} is set to \code{TRUE}, bootstrapping is carried out. This is done by resampling, with replacement, the rows of \code{data} and then computing the alpha statistic for the resulting matrix. The elements of argument \code{control} are used to control the bootstrap computation.
#'
#' @param data a matrix of scores. Each row corresponds to a unit, each column a coder.
#' @param level the level of measurement, one of \code{"nominal"}, \code{"ordinal"}, \code{"interval"}, or \code{"ratio"}; or a user-defined distance function.
#' @param confint logical; if \code{TRUE}, a bootstrap sample is produced.
#' @param verbose logical; if \code{TRUE}, various messages are printed to the console. Note that if \code{confint = TRUE} a progress bar (\code{\link[pbapply]{pblapply}}) is displayed (if possible) during the bootstrap computation.
#' @param control a list of control parameters.
#'    \describe{
#'        \item{\code{bootit}}{the size of the bootstrap sample. This applies when \code{confint = TRUE}. Defaults to 1,000.}
#'        \item{\code{nodes}}{the desired number of nodes in the cluster.}
#'        \item{\code{parallel}}{logical; if \code{TRUE} (the default), bootstrapping is done in parallel.}
#'        \item{\code{type}}{one of the supported cluster types for \code{\link[parallel]{makeCluster}}. Defaults to \code{"SOCK"}.}
#' }
#'
#' @return Function \code{krippendorffs.alpha} returns an object of class \code{"krippendorffsalpha"}, which is a list comprising the following elements.
#'         \item{boot.sample}{when applicable, the bootstrap sample.}
#'         \item{call}{the matched call.}
#'         \item{coders}{the number of coders.}
#'         \item{alpha.hat}{the estimate of alpha.}
#'         \item{confint}{the value of argument \code{confint}.}
#'         \item{control}{the list of control parameters.}
#'         \item{data}{the matrix of scores.}
#'         \item{D.e}{the estimate of total variation.}
#'         \item{D.o}{the estimate of within-unit variation.}
#'         \item{level}{the level of measurement.}
#'         \item{units}{the number of units.}
#'         \item{verbose}{the value of argument \code{verbose}.}
#'
#' @references
#' Krippendorff, K. (2013). Computing Krippendorff's alpha-reliability. Technical report, University of Pennsylvania.
#'
#' @export
#'
#' @examples
#' # The following data were presented in Krippendorff (2013).
#'
#' nominal = matrix(c(1,2,3,3,2,1,4,1,2,NA,NA,NA,
#'                    1,2,3,3,2,2,4,1,2,5,NA,3,
#'                    NA,3,3,3,2,3,4,2,2,5,1,NA,
#'                    1,2,3,3,2,4,4,1,2,5,1,NA), 12, 4)
#' nominal
#' fit.nom = krippendorffs.alpha(nominal, level = "nominal", confint = TRUE, verbose = TRUE,
#'                               control = list(bootit = 100, parallel = FALSE))
#' summary(fit.nom)
#' confint(fit.nom, level = 0.99)

krippendorffs.alpha = function(data, level = c("interval", "nominal", "ordinal", "ratio"), confint = TRUE,
                               verbose = FALSE, control = list())
{
    call = match.call()
    if (missing(data) || ! is.matrix(data) || ! is.numeric(data))
        stop("You must supply a numeric data matrix.")
    if (! is.character(level) && ! is.function(level))
        stop("'level' must be a distance function or one of \"interval\", \"nominal\", \"ordinal\", or \"ratio\".")
    else if (is.function(level))
    {
        dist = level
        level = "user specified distance function"
    }
    else
    {
        if (level == "interval" || level == "ordinal")
            dist = interval.dist
        else if (level == "nominal")
            dist = nominal.dist
        else if (level == "ratio")
            dist = ratio.dist
        else
            stop("'level' must be a distance function or one of \"interval\", \"nominal\", \"ordinal\", or \"ratio\".")
    }
    if (! is.logical(confint) || length(confint) > 1)
        stop("'confint' must be a logical value.")
    if (! is.logical(verbose) || length(verbose) > 1)
        stop("'verbose' must be a logical value.")
    if (! is.list(control))
        stop("'control' must be a list.")
    control = krippendorff.control(confint, verbose, control)
    D.o = krippendorff.within(data, dist)
    D.e = krippendorff.total(data, dist)
    alpha.hat = 1 - D.o / D.e
    names(alpha.hat) = "alpha"
    object = list(units = nrow(data), coders = ncol(data), data = data, level = level, confint = confint,
                  verbose = verbose, alpha.hat = alpha.hat, control = control, call = call, D.o = D.o,
                  D.e = D.e, dist = dist)
    if (confint)
    {
        boot.sample = krippendorff.bootstrap(object)
        object$boot.sample = boot.sample
    }
    class(object) = "krippendorffsalpha"
    object
}

#' Compute the squared difference between two scores.
#'
#' @details This function computes the squared difference between two scores. This may be an appropriate distance function for the interval level of measurement. \code{NA}'s are handled gracefully.
#' @param x a score.
#' @param y a score.
#' @return \eqn{(x-y)^2}, or 0 if \code{x} or \code{y} is \code{NA}.
#' @seealso \code{\link{nominal.dist}}, \code{\link{ratio.dist}}
#' @export
    
interval.dist = function(x, y)
{
    d = (x - y)^2
    if (is.na(d))
        d = 0
    d
}

#' Apply the discrete metric to two scores.
#'
#' @details This function applies the discrete metric to two scores. This may be an appropriate distance function for the nominal level of measurement. \code{NA}'s are handled gracefully.
#' @param x a score.
#' @param y a score.
#' @return 0 if \code{x} is equal to \code{y} or if either is \code{NA}, 1 otherwise.
#' @seealso \code{\link{interval.dist}}, \code{\link{ratio.dist}}
#' @export

nominal.dist = function(x, y)
{
    temp = (x == y)
    if (is.na(temp) || temp)
        d = 0
    else
        d = 1
    d
}

#' Apply a ratio distance function to two scores.
#'
#' @details This function applies a ratio distance function to two scores. This may be an appropriate distance function for the ratio level of measurement. \code{NA}'s are handled gracefully.
#' @param x a score.
#' @param y a score.
#' @return \eqn{(x-y)^2/(x+y)^2}, or 0 if \code{x} or \code{y} is \code{NA}.
#' @seealso \code{\link{interval.dist}}, \code{\link{nominal.dist}}
#' @export

ratio.dist = function(x, y)
{
    d = (x - y)^2 / (x + y)^2
    if (is.na(d))
        d = 0
    d
}

#' Compute DFBETAs for units and/or coders.
#'
#' @details This function computes DFBETAS for one or more units and/or one or more coders.
#'
#' @param model a fitted model object, the result of a call to \code{\link{krippendorffs.alpha}}.
#' @param units a vector of integers. A DFBETA will be computed for each of the corresponding units.
#' @param coders a vector of integers. A DFBETA will be computed for each of the corresponding coders.
#' @param ... additional arguments. These are ignored.
#
#' @return A list comprising at most two elements.
#'         \item{dfbeta.units}{a vector containing DFBETAS for the units specified via argument \code{units}.}
#'         \item{dfbeta.coders}{a vector containing DFBETAS for the coders specified via argument \code{coders}.}
#'
#' @method influence krippendorffsalpha
#'
#' @references
#' Young, D. S. (2017). \emph{Handbook of Regression Methods}. CRC Press.
#' @references
#' Krippendorff, K. (2013). Computing Krippendorff's alpha-reliability. Technical report, University of Pennsylvania.
#'
#' @export
#'
#' @examples
#' # The following data were presented in Krippendorff (2013).
#'
#' nominal = matrix(c(1,2,3,3,2,1,4,1,2,NA,NA,NA,
#'                    1,2,3,3,2,2,4,1,2,5,NA,3,
#'                    NA,3,3,3,2,3,4,2,2,5,1,NA,
#'                    1,2,3,3,2,4,4,1,2,5,1,NA), 12, 4)
#' fit.nom = krippendorffs.alpha(nominal, level = "nominal", confint = FALSE)
#' summary(fit.nom)
#' (inf = influence(fit.nom, units = c(6, 11), coders = c(2, 3)))

influence.krippendorffsalpha = function(model, units, coders, ...)
{
    dfbeta.units = NULL
    if (! missing(units))
    {
        dfbeta.units = numeric(length(units))
        k = 1
        for (j in units)
        {
            fit = krippendorffs.alpha(model$data[-j, ], level = model$level, confint = FALSE, control = model$control)
            dfbeta.units[k] = model$alpha.hat - fit$alpha.hat
            k = k + 1
        }
        names(dfbeta.units) = units
    }
    dfbeta.coders = NULL
    if (! missing(coders))
    {
        dfbeta.coders = numeric(length(coders))
        k = 1
        for (j in coders)
        {
            fit = krippendorffs.alpha(model$data[, -j], level = model$level, confint = FALSE, control = model$control)
            dfbeta.coders[k] = model$alpha.hat - fit$alpha.hat
            k = k + 1
        }
        names(dfbeta.coders) = coders
    }
    result = list()
    if (! is.null(dfbeta.units))
        result$dfbeta.units = dfbeta.units
    if (! is.null(dfbeta.coders))
        result$dfbeta.coders = dfbeta.coders
    result        
}

#' Print a summary of a Krippendorff's Alpha fit.
#'
#' @details This function prints a summary of the fit. First the values of the control parameters (defaults and/or values supplied in the call) are printed. Then a table of estimates is shown. If applicable, the table includes confidence intervals.
#'
#' @param object an object of class \code{"krippendorffsalpha"}, the result of a call to \code{\link{krippendorffs.alpha}}.
#' @param conf.level the confidence level for the confidence intervals. The default is 0.95.
#' @param digits the number of significant digits to display. The default is 4.
#' @param \dots additional arguments. These are passed to \code{\link{quantile}}.
#'
#' @seealso \code{\link{krippendorffs.alpha}}
#'
#' @method summary krippendorffsalpha
#'
#' @references
#' Nissi, M. J., Mortazavi, S., Hughes, J., Morgan, P., and Ellermann, J. (2015). T2* relaxation time of acetabular and femoral cartilage with and without intra-articular Gd-DTPA2 in patients with femoroacetabular impingement. \emph{American Journal of Roentgenology}, \bold{204}(6), W695.
#'
#' @export
#'
#' @examples
#' # Fit a subset of the cartilage data. Compute bootstrap confidence intervals
#' # using a bootstrap sample size of 1,000. Display a summary of the results,
#' # including a 99% confidence interval. Also plot the results.
#'
#' data(cartilage)
#' cartilage = as.matrix(cartilage[1:100, ])
#' fit.cart = krippendorffs.alpha(cartilage, level = "interval", confint = TRUE,
#'                                control = list(bootit = 1000, parallel = FALSE))
#' summary(fit.cart, conf.level = 0.99)
#' dev.new()
#' plot(fit.cart, xlim = c(0.7, 0.9), xlab = "Bootstrap Estimates",
#'      main = "Results for Cartilage Data")

summary.krippendorffsalpha = function(object, conf.level = 0.95, digits = 4, ...)
{
	cat("\nKrippendorff's Alpha\n\n")
	cat("Data:", object$units, "units x", object$coders, "coders\n")
	ans = NULL
	ans$units = object$units
	ans$coders = object$coders
    cat("\nCall:\n\n")
    print(object$call)
	ans$call = object$call
    cat("\nControl parameters:\n")
    if (length(object$control) > 0)
    {
        control.table = cbind(unlist(c(object$control, "")))
        colnames(control.table) = ""
		ans$control.table = control.table
        print(control.table, quote = FALSE)
    }
    else
        cat("\nNA\n\n")
    if (! object$confint)
    {
        confint = matrix(rep(NA, 2), ncol = 2)
        coef.table = cbind(object$alpha.hat, confint)
    }
    else
    {
        boot.sample = object$boot.sample
        if (! is.null(boot.sample))
        {
            lo = (1 - conf.level) / 2
            probs = c(lo, 1 - lo)
            confint = quantile(boot.sample, probs, ...)
        }
        coef.table = cbind(object$alpha.hat, t(confint))
    }
    colnames(coef.table) = c("Estimate", "Lower", "Upper")
    rownames(coef.table) = names(object$alpha.hat)
	ans$coef.table = coef.table
    cat("Results:\n\n")
    print(signif(coef.table, digits = digits))
    cat("\n")
	invisible(ans)
}

#' Compute a confidence interval for Krippendorff's Alpha.
#'
#' @details This function computes a bootstrap confidence interval for alpha, assuming that \code{\link{krippendorffs.alpha}} was called with \code{confint = TRUE}.
#'
#' @param object an object of class \code{"krippendorffsalpha"}, the result of a call to \code{\link{krippendorffs.alpha}}.
#' @param parm always ignored since there is only one parameter.
#' @param level the desired confidence level for the interval. The default is 0.95.
#' @param \dots additional arguments. These are passed to \code{\link{quantile}}.
#' @return A vector with entries giving lower and upper confidence limits. These will be labelled as (1-level)/2 and 1 - (1-level)/2.
#' @seealso \code{\link{krippendorffs.alpha}}
#'
#' @method confint krippendorffsalpha
#'
#' @references
#' Nissi, M. J., Mortazavi, S., Hughes, J., Morgan, P., and Ellermann, J. (2015). T2* relaxation time of acetabular and femoral cartilage with and without intra-articular Gd-DTPA2 in patients with femoroacetabular impingement. \emph{American Journal of Roentgenology}, \bold{204}(6), W695.
#'
#' @export
#'
#' @examples
#' # Fit a subset of the cartilage data. Compute bootstrap confidence intervals
#' # using a bootstrap sample size of 1,000.
#'
#' data(cartilage)
#' cartilage = as.matrix(cartilage[1:100, ])
#' fit.cart = krippendorffs.alpha(cartilage, level = "interval", confint = TRUE,
#'                                control = list(bootit = 1000, parallel = FALSE))
#' fit.cart$alpha.hat
#' confint(fit.cart, level = 0.99)

confint.krippendorffsalpha = function(object, parm = "alpha", level = 0.95, ...)
{
    if (! object$confint)
        return(rep(NA, 2))
    else
    {
        boot.sample = object$boot.sample
        if (! is.null(boot.sample))
        {
            lo = (1 - level) / 2
            hi = 1 - lo
            probs = c(lo, hi)
            confint = quantile(boot.sample, probs, ...)
            names(confint) = probs
        }
    }
    confint
}

#' Plot the results of a Krippendorff's Alpha analysis.
#'
#' @details This function plots the results of a Krippendorff's Alpha analysis, assuming that \code{\link{krippendorffs.alpha}} was called with \code{confint = TRUE}.
#'
#' @param x an object of class \code{"krippendorffsalpha"}, the result of a call to \code{\link{krippendorffs.alpha}}.
#' @param y always ignored.
#' @param level the desired confidence level for the interval. The default is 0.95.
#' @param type the method used to compute sample quantiles. This argument is passed to \code{\link{quantile}}. The default is 7.
#' @param density logical; if \code{TRUE}, a kernel density estimate is plotted.
#' @param lty.density the line type for the kernel density estimate. The default is 1.
#' @param lty.estimate the line type for the estimate of alpha. The default is 1.
#' @param lty.interval the line type for the confidence limits. The default is 2.
#' @param col.density the color for the kernel density estimate. The default is black.
#' @param col.estimate the color for the estimate of alpha. The default is orange.
#' @param col.interval the color for the confidence limits. The default is blue.
#' @param lwd.density the line width for the kernel density estimate. The default is 3.
#' @param lwd.estimate the line width for the estimate of alpha. The default is 3.
#' @param lwd.interval the line width for the confidence limits. The default is 3.
#' @param \dots additional arguments. These are passed to \code{\link{hist}}.
#'
#' @method plot krippendorffsalpha
#'
#' @references
#' Krippendorff, K. (2013). Computing Krippendorff's alpha-reliability. Technical report, University of Pennsylvania.
#'
#' @export
#'
#' @examples
#' # The following data were presented in Krippendorff (2013).
#'
#' nominal = matrix(c(1,2,3,3,2,1,4,1,2,NA,NA,NA,
#'                    1,2,3,3,2,2,4,1,2,5,NA,3,
#'                    NA,3,3,3,2,3,4,2,2,5,1,NA,
#'                    1,2,3,3,2,4,4,1,2,5,1,NA), 12, 4)
#' fit.nom = krippendorffs.alpha(nominal, level = "nominal", confint = TRUE, verbose = TRUE,
#'                               control = list(bootit = 1000, parallel = FALSE))
#' dev.new()
#' plot(fit.nom, main = "Results for Nominal Data", xlab = "Bootstrap Estimates", density = FALSE)

plot.krippendorffsalpha = function(x, y = NULL, level = 0.95, type = 7, density = TRUE,
	                               lty.density = 1, lty.estimate = 1, lty.interval = 2,
	                               col.density = "black", col.estimate = "orange", col.interval = "blue",
								   lwd.density = 3, lwd.estimate = 3, lwd.interval = 3, ...)
{
    if (! x$confint)
        stop("There is nothing to plot. Call function krippendorffs.alpha with confint = TRUE.")
    else
    {
        boot.sample = x$boot.sample
        if (! is.null(boot.sample))
        {
            lo = (1 - level) / 2
            hi = 1 - lo
            probs = c(lo, hi)
            confint = quantile(boot.sample, probs, type = type)
        }
		hist(boot.sample, prob = TRUE, ...)
		if (density)
		    lines(density(boot.sample), lty = lty.density, col = col.density, lwd = lwd.density)
		abline(v = x$alpha.hat, lty = lty.estimate, col = col.estimate, lwd = lwd.estimate)
		abline(v = confint, lty = lty.interval, col = col.interval, lwd = lwd.interval)
    }
}
