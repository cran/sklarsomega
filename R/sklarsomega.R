
#' Build a Sklar's Omega correlation matrix.
#'
#' @details This function accepts a data matrix and uses its column names to build the appropriate block-diagonal correlation matrix. If gold-standard scores are included, they should be in the first column of the matrix, and that column should be named 'g'. For a given coder, the column name should begin with 'c', and then the coder number and score number should follow, separated by '.' (so that multi-digit numbers can be accommodated). For example, 'c.12.2' denotes the second score for coder 12.
#'
#' See the package vignette for detailed information regarding the structure of the correlation matrix.
#'
#' Note that this function is called by \code{\link{sklars.omega}} and so is not a user-level
#' function, per se. We expose the function so that interested users can more easily carry out
#' simulation studies.
#'
#' @param data a matrix of scores. Each row corresponds to a unit, each column a coder.
#'
#' @return \code{build.R} returns a list comprising two elements.
#'         \item{R}{the correlation matrix.}
#'         \item{onames}{a character vector that contains names for the parameters of the correlation matrix.}
#'
#' @seealso \code{\link{sklars.omega}}, \code{\link{check.colnames}}
#'
#' @export

build.R = function(data)
{
    cnames = colnames(data)
    cols = ncol(data)
    params = 1
    coders = hash()
    onames = NULL
    if (cnames[1] == "g")
    {
        params = params + 1
        coders[["g"]] = 1
        onames = c(onames, "gold")
    }
    else
    {
        current = strsplit(cnames[1], ".", fixed = TRUE)[[1]][-1]
        coders[["1"]] = 1
    }
    for (j in 2:cols)
    {
        current = strsplit(cnames[j], ".", fixed = TRUE)[[1]][-1]
        if (is.null(coders[[current[1]]]))
            coders[[current[1]]] = 1
        else
            coders[[current[1]]] = coders[[current[1]]] + 1
    }
    vals = values(coders)
    for (j in 1:length(vals))
        if (vals[j] > 1)
            params = params + 1
    omega = seq(0.1, 0.9, length = params)
    block = diag(1, cols)
    if (! is.null(coders[["g"]]))
    {
        block[1, 2:cols] = block[2:cols, 1] = omega[1]
        block[block == 0] = omega[2]
        if (params > 2 || length(vals) > 2)
            onames = c(onames, "inter")
        o = 3
        k = 2
    }
    else
    {
        block[block != 1] = omega[1]
        o = 2
        k = 1
        if (length(vals) > 1)
            onames = c(onames, "inter")
    }
    if (length(vals) > 1)
    {
        for (j in 1:length(vals))
        {
            if (vals[j] > 1)
            {
                sub = matrix(omega[o], vals[j], vals[j])
                diag(sub) = 1
                block[k:(k + vals[j] - 1), k:(k + vals[j] - 1)] = sub
                k = k + vals[j] - 1
                o = o + 1
                onames = c(onames, "intra")
            }
            else
                k = k + 1
        }
    }
    else
        onames = c(onames, "intra")
    block.sizes = apply(data, 1, function(row) { sum(! is.na(row)) })
    dat = data[block.sizes > 1, ]
    block.sizes = block.sizes[block.sizes > 1]
    block.list = NULL
    for (i in 1:nrow(dat))
    {
        present = which(! is.na(dat[i, ]))
        block.list[[i]] = block[present, present]
    }
    R = as.matrix(bdiag(block.list))
    list(R = R, onames = onames)
}


objective.ML = function(theta, y, R, dist = c("gaussian", "laplace", "t", "beta", "gamma", "empirical"), F.hat = NULL)
{
    vals = sort(unique(R[R != 0 & R != 1]))
    m = length(vals)
    omega = theta[1:m]
    for (j in 1:m)
        R[R == vals[j]] = omega[j]
    rest = theta[(m + 1):length(theta)]
    R = as.spam(R)
    Rinv = try(solve.spam(R), silent = TRUE)
    if (inherits(Rinv, "try-error"))
        return(1e6)
    dist = match.arg(dist)
    u = switch(dist,
               gaussian = pnorm(y, mean = rest[1], sd = rest[2]),
               laplace = plaplace(y, location = rest[1], scale = rest[2]),
               t = pt(y, df = rest[1], ncp = rest[2]),
               beta = pbeta(y, shape1 = rest[1], shape2 = rest[2]),
               gamma = pgamma(y, shape = rest[1], rate = rest[2]),
               empirical = F.hat(y))
    z = qnorm(u)
    z = ifelse(z == -Inf, qnorm(0.0001), z)
    z = ifelse(z == Inf, qnorm(0.9999), z)
    qform =  try(t(z) %*% Rinv  %*% z, silent = TRUE)
    if (inherits(qform, "try-error"))
        return(1e6)
    if (dist == "empirical")
        z = 0
    logf = switch(dist,
                  gaussian = dnorm(y, mean = rest[1], sd = rest[2], log = TRUE),
                  laplace = dlaplace(y, location = rest[1], scale = rest[2], log = TRUE),
                  t = dt(y, df = rest[1], ncp = rest[2], log = TRUE),
                  beta = dbeta(y, shape1 = rest[1], shape2 = rest[2], log = TRUE),
                  gamma = dgamma(y, shape = rest[1], rate = rest[2], log = TRUE),
                  empirical = 0)
    obj = as.numeric(0.5 * (determinant(R, logarithm = TRUE)$modulus + qform - sum(z^2)) - sum(logf))
    obj
}


#' Compute the cumulative distribution function for a categorical distribution.
#'
#' @details This function uses \code{\link[LaplacesDemon]{dcat}} to compute the cdf for the categorical distribution with support \eqn{1, \dots , K} and probabilities \eqn{p = (p_1, \dots , p_K)'}.
#'
#' @param q vector of quantiles.
#' @param p vector of probabilities. 
#
#' @return \code{pcat} returns a probability or a vector of probabilities, depending on the dimension of \code{q}.
#'
#' @export

pcat = function(q, p)
{
    n = length(q)
    pr = numeric(n)
    for (j in 1:n)
        if (q[j] > 0)
            pr[j] = sum(LaplacesDemon::dcat(1:q[j], p = p))
    pr
}


#' Compute the quantile function for a categorical distribution.
#'
#' @details This function computes quantiles for the categorical distribution with support \eqn{1, \dots , K} and probabilities \eqn{p = (p_1, \dots , p_K)'}.
#'
#' @param pr vector of probabilities.
#' @param p vector of proabilities.
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X \le  x)}, otherwise, \eqn{P(X > x)}.
#' @param log.pr logical; if TRUE, probabilities \eqn{pr} are handled as log\eqn{(pr)}.
#
#' @return \code{qcat} returns a vector of quantiles.
#'
#' @export

qcat = function (pr, p, lower.tail = TRUE, log.pr = FALSE) 
{
    if (! is.vector(pr)) 
        pr = as.vector(pr)
    if (! is.vector(p)) 
        p = as.vector(p)
    p = p / sum(p)
    if (log.pr == FALSE)
    {
        if (any(pr < 0) | any(pr > 1)) 
            stop("pr must be in [0,1].")
    }
    else if (any(! is.finite(pr)) | any(pr > 0)) 
        stop("pr, as a log, must be in (-Inf,0].")
    if (sum(p) != 1) 
        stop("sum(p) must be 1.")
    if (lower.tail == FALSE) 
        pr = 1 - pr
    breaks = c(0, cumsum(p))
    if (log.pr == TRUE) 
        breaks = log(breaks)
    breaks = matrix(breaks, length(pr), length(breaks), byrow = TRUE)
    x = rowSums(pr > breaks)
    x
}


objective.DT = function(theta, y, R)
{
    vals = sort(unique(R[R != 0 & R != 1]))
    m = length(vals)
    omega = theta[1:m]
    for (j in 1:m)
        R[R == vals[j]] = omega[j]
    p = theta[(m + 1):length(theta)]
    p = p / sum(p)
    R = as.spam(R)
    Rinv = try(solve.spam(R), silent = TRUE)
    if (inherits(Rinv, "try-error"))
        return(1e6)
    u = (sklarsomega::pcat(y, p) + sklarsomega::pcat(y - 1, p)) / 2
    z = qnorm(u)
    z = ifelse(z == -Inf, qnorm(0.0001), z)
    z = ifelse(z == Inf, qnorm(0.9999), z)
    qform =  try(t(z) %*% Rinv  %*% z, silent = TRUE)
    if (inherits(qform, "try-error"))
        return(1e6)
    logf = dcat(y, p, log = TRUE)
    as.numeric(0.5 * (determinant(R, logarithm = TRUE)$modulus + qform - sum(z^2)) - sum(logf))
}


pbivnormf = function (K, omega = 0)
{
    correl = rep(omega, nrow(K))
    lower = as.double(c(0, 0))
    infin = as.integer(c(0, 0))
    uppera = as.double(K[, 1])
    upperb = as.double(K[, 2])
    lt = as.integer(nrow(K))
    prob = double(lt)
    correl = as.double(correl)
    ans = .Fortran("PBIVNORM", prob, lower, uppera, upperb, infin, correl, lt, PACKAGE = "sklarsomega")[[1]]
    ans
}


neighbor.list = function(R) 
{
    n = nrow(R)
    N = vector("list", n)
    for (i in 1:n)
    {
        temp = which(as.logical(R[i, ]))
        N[[i]] = temp[temp > i]
    }
    N
}


objective.CML = function(theta, y, R, N, K)
{
    vals = sort(unique(R[R != 0 & R != 1]))
    m = length(vals)
    omega = theta[1:m]
    for (j in 1:m)
        R[R == vals[j]] = omega[j]
    p = theta[(m + 1):length(theta)]
    p = p / sum(p)
    u.0 = sklarsomega::pcat(y, p)
    u.1 = sklarsomega::pcat(y - 1, p)
    u.0 = ifelse(u.0 < 0, 0, u.0)
    u.1 = ifelse(u.1 < 0, 0, u.1)
    u.0 = ifelse(u.0 > 1, 1, u.0)
    u.1 = ifelse(u.1 > 1, 1, u.1)
    z.0 = qnorm(u.0)
    z.1 = qnorm(u.1)
    z.0 = ifelse(z.0 == -Inf, qnorm(0.0001), z.0)
    z.1 = ifelse(z.1 == -Inf, qnorm(0.0001), z.1)
    z.0 = ifelse(z.0 == Inf, qnorm(0.9999), z.0)
    z.1 = ifelse(z.1 == Inf, qnorm(0.9999), z.1)
    obj = 0
    u = c(1, -1, -1, 1)
    for (i in 1:length(N))
    {
        for (j in N[[i]])
        {
            K[, 1] = c(z.0[i], z.0[i], z.1[i], z.1[i])
            K[, 2] = c(z.0[j], z.1[j], z.0[j], z.1[j])
            s = pbivnormf(K, omega = R[i, j])
            s = t(s) %*% u
            if (! is.na(s) && s > 0)
                obj = obj + log(s)
        }
    }
    -obj
}


is.positiveinteger = function(x, tol = .Machine$double.eps^0.5)
{
    x = suppressWarnings(as.numeric(x))
    if (is.na(x))
        return(FALSE)
    (x > 0 & abs(x - round(x)) < tol)
}


#' Check the column names of a Sklar's Omega data matrix for correctness.
#'
#' @details This function performs a somewhat rudimentary validation of the column names. At most one column may be labeled 'g', and said column must be the first. The only other valid format is 'c.C.S', where 'C' denotes coder number and 'S' denotes the Sth score for coder C (both positive whole numbers). It is up to the user to ensure that the coder and score indices make sense and are ordered correctly, i.e., coders, and scores for a given coder, are numbered consecutively.
#'
#' @param data a matrix of scores. Each row corresponds to a unit, each column a coder.
#
#' @return \code{check.colnames} returns a list comprising two elements.
#'         \item{success}{logical; if TRUE, the column names passed the test.}
#'         \item{cols}{if \code{success} is FALSE, vector \code{cols} contains the numbers of the problematic column names.}
#'
#' @references
#' Krippendorff, K. (2013). Computing Krippendorff's alpha-reliability. Technical report, University of Pennsylvania.
#'
#' @export
#'
#' @examples
#' # The following data were presented in Krippendorff (2013).
#'
#' data = matrix(c(1,2,3,3,2,1,4,1,2,NA,NA,NA,
#'                 1,2,3,3,2,2,4,1,2,5,NA,3,
#'                 NA,3,3,3,2,3,4,2,2,5,1,NA,
#'                 1,2,3,3,2,4,4,1,2,5,1,NA), 12, 4)
#' colnames(data) = c("c.1.1", "c.2.1", "c.3.1", "c.4.1")
#' data
#' (check.colnames(data))
#'
#' # Introduce errors for columns 1 and 4.
#'
#' colnames(data) = c("c.a.1", "c.2.1", "c.3.1", "C.4.1")
#' (check.colnames(data))
#'
#' # The following scenario passes the check but is illogical.
#'
#' colnames(data) = c("g", "c.2.1", "c.1.47", "c.2.1")
#' (check.colnames(data))

check.colnames = function(data)
{
    cnames = colnames(data)
    m = ncol(data)
    success = TRUE
    cols = NULL
    temp = strsplit(cnames, ".", fixed = TRUE)
    j = 1
    if (length(temp[[1]]) == 1)
    {
        if (temp[[1]] != "g")
        {
            success = FALSE
            cols = c(cols, 1)
        }
        j = 2
    }
    while (j <= m)
    {
        current = temp[[j]]
        if (length(current) != 3 |
            current[1] != "c" |
            ! is.positiveinteger(current[2]) |
            ! is.positiveinteger(current[3]))
        {
            success = FALSE
            cols = c(cols, j)
        }
        j = j + 1
    }
    obj = list(success = success, cols = cols)
}


sklarsomega.control = function(level, confint, verbose, control)
{
    if (length(control) != 0)
    {
        nms = match.arg(names(control), c("bootit", "parallel", "type", "nodes", "dist"), several.ok = TRUE)
        control = control[nms]
    }
    boot = FALSE
    if (level == "interval")
    {
        dist = control$dist
        if (is.null(dist) || length(dist) > 1 || ! dist %in% c("gaussian", "laplace", "t", "gamma", "empirical"))
            stop("\nControl parameter 'dist' must be \"gaussian\", \"laplace\", \"t\", \"gamma\", or \"empirical\".")
        if (confint == "bootstrap")
            boot = TRUE
    }
    else if (level == "ratio")
    {
        control$dist = "beta"
        if (confint == "bootstrap")
            boot = TRUE
    }
    else if (level %in% c("nominal", "ordinal"))
    {
        control$dist = "categorical"
        if (confint != "none")
            boot = TRUE
    }
    if (boot)
    {
        bootit = control$bootit
        if (is.null(bootit) || ! is.numeric(bootit) || length(bootit) > 1 || bootit != as.integer(bootit) || bootit < 2)
        {
            if (verbose)
                cat("\nControl parameter 'bootit' must be a positive integer > 1. Setting it to the default value of 1,000.\n")
            control$bootit = 1000
        }
        if (is.null(control$parallel) || ! is.logical(control$parallel) || length(control$parallel) > 1)
        {
            if (verbose)
                cat("\nControl parameter 'parallel' must be a logical value. Setting it to the default value of FALSE.\n")
            control$parallel = FALSE
            control$type = control$nodes = NULL
        }
        if (control$parallel)
        {
            if (requireNamespace("parallel", quietly = TRUE))
            {
                if (is.null(control$type) || length(control$type) > 1 || ! control$type %in% c("SOCK", "PVM", "MPI", "NWS"))
                {
                    if (verbose)
                        cat("\nControl parameter 'type' must be \"SOCK\", \"PVM\", \"MPI\", or \"NWS\". Setting it to \"SOCK\".\n")
                    control$type = "SOCK"
                }
                nodes = control$nodes
                if (is.null(control$nodes) || ! is.numeric(nodes) || length(nodes) > 1 || nodes != as.integer(nodes) || nodes < 2)
                    stop("Control parameter 'nodes' must be a whole number greater than 1.")
            }
            else
            {
                if (verbose)
                    cat("\nParallel computation requires package parallel. Setting control parameter 'parallel' to FALSE.\n")
                control$parallel = FALSE
                control$type = control$nodes = NULL 
            }
        }
    }
    else
        control$bootit = control$parallel = control$type = control$nodes = NULL
    control
}


bootstrap.helper = function(object)
{
    boot.sample = matrix(NA, object$control$bootit, object$npar)
    sim = simulate(object, object$control$bootit)
    data = vector("list", object$control$bootit)
    for (j in 1:object$control$bootit)
    {
        y = sim[, j]
        dat = t(object$data)
        dat[! is.na(dat)] = y
        dat = t(dat)
        data[[j]] = dat
    }
    rm(sim)
    if (! object$control$parallel)
    {
        if (object$verbose && requireNamespace("pbapply", quietly = TRUE))
        {
            cat("\n")
            gathered = pbapply::pblapply(data, sklars.omega, level = object$level, confint = "none", control = object$control)
            for (j in 1:object$control$bootit)
            {
                fit = gathered[[j]]
                if (inherits(fit, "try-error") || fit$convergence != 0)
                    warning(fit$message)
                else
                    boot.sample[j, ] = fit$coef
            }
        }
        else
        {
            for (j in 1:object$control$bootit)
            {
                fit = try(suppressWarnings(sklars.omega(data[[j]], level = object$level, confint = "none", control = object$control)), silent = TRUE)
                if (inherits(fit, "try-error") || fit$convergence != 0)
                    warning(fit$message)
                else
                    boot.sample[j, ] = fit$coef
            }
        }
    }
    else
    {
        cl = parallel::makeCluster(object$control$nodes, object$control$type)
        parallel::clusterEvalQ(cl, library(sklarsomega))
        if (object$verbose && requireNamespace("pbapply", quietly = TRUE))
        {
            cat("\n")
            gathered = pbapply::pblapply(data, sklars.omega, level = object$level, confint = "none", control = object$control, cl = cl)
        }
        else
            gathered = parallel::clusterApplyLB(cl, data, sklars.omega, level = object$level, confint = "none", control = object$control)
        parallel::stopCluster(cl)
        for (j in 1:object$control$bootit)
        {
            fit = gathered[[j]]
            if (inherits(fit, "try-error") || fit$convergence != 0)
                warning(fit$message)
            else
                boot.sample[j, ] = fit$coef
        }
    }
    boot.sample = boot.sample[complete.cases(boot.sample), ]
    if (is.vector(boot.sample))
        boot.sample = as.matrix(boot.sample)
    boot.sample
}


grad.helper = function(y, params, R, N, K)
{
    if (is.null(N))
        gr = try(-grad(objective.DT, params, y = y, R = R), silent = TRUE)
    else
        gr = try(-grad(objective.CML, params, y = y, R = R, N, K), silent = TRUE)
    gr
}


sandwich.helper = function(object)
{
    N = K = NULL
    if (object$method == "CML")
    {
        N = neighbor.list(object$R)
        K = matrix(NA, 4, 2)
    }
    sim = simulate(object, object$control$bootit)
    data = vector("list", object$control$bootit)
    for (j in 1:object$control$bootit)
    {
        y = sim[, j]
        data[[j]] = y
    }
    rm(sim)
    meat = 0
    if (! object$control$parallel)
    {
        if (object$verbose && requireNamespace("pbapply", quietly = TRUE))
        {
            cat("\n")
            gathered = pbapply::pblapply(data, grad.helper, params = object$coef, R = object$R, N = N, K = K)
            b = 1
            for (j in 1:object$control$bootit)
            {
                gr = gathered[[j]]
                if (! inherits(gr, "try-error"))
                {
                    meat = meat + gr %o% gr
                    b = b + 1
                }
            }
            meat = meat / b
        }
        else
        {
            b = 1
            for (j in 1:object$control$bootit)
            {
                gr = try(grad.helper(data[[j]], params = object$coef, R = object$R, N = N, K = K), silent = TRUE)
                if (! inherits(gr, "try-error"))
                {
                    meat = meat + gr %o% gr
                    b = b + 1
                }
            }
            meat = meat / b
        }
    }
    else
    {
        cl = parallel::makeCluster(object$control$nodes, object$control$type)
        parallel::clusterEvalQ(cl, library(sklarsomega))
        if (object$verbose && requireNamespace("pbapply", quietly = TRUE))
        {
            cat("\n")
            gathered = pbapply::pblapply(data, grad.helper, params = object$coef, R = object$R, N = N, K = K, cl = cl)
        }
        else
            gathered = parallel::clusterApplyLB(cl, data, grad.helper, params = object$coef, R = object$R, N = N, K = K)
        parallel::stopCluster(cl)
        b = 1
        for (j in 1:object$control$bootit)
        {
            gr = gathered[[j]]
            if (! inherits(gr, "try-error"))
            {
                meat = meat + gr %o% gr
                b = b + 1
            }
        }
        meat = meat / b
    }
    meat
}


#' Apply Sklar's Omega.
#'
#' @details This is the package's flagship function. It applies the Sklar's Omega methodology to nominal, ordinal, interval, or ratio outcomes, and, if desired, produces confidence intervals. Parallel computing is supported, when applicable, and other measures (e.g., sparse matrix operations) are taken in the interest of computational efficiency.
#'
#' If the level of measurement is nominal or ordinal, the scores (which must take values in \eqn{1, \dots , K}) are assumed to share a common categorical marginal distribution. If \eqn{K} is less than 5, a composite marginal likelihood (CML) approach is used. If \eqn{K} is greater than or equal to 5, the distributional transform (DT) approximation is used. In either case, two types of confidence interval are available: bootstrap or asymptotic. See the package vignette for details.
#'
#' If the level of measurement is ratio, the scores are assumed to share a common beta distribution, and the method of maximum likelihood (ML) is used. Both bootstrap and asymptotic confidence intervals are available.
#'
#' If the level of measurement is interval, control parameter \code{dist} must be used to select a marginal distribution from among \code{"gaussian"}, \code{"laplace"}, \code{"t"}, \code{"gamma"}, and \code{"empirical"}. The ML method is used unless \code{dist = "empirical"}, in which case conditional maximum likelihood is used, i.e., the copula parameters are estimated conditional on the sample distribution function of the scores. For the ML method, both bootstrap and asymptotic confidence intervals are available. When \code{dist = "empirical"}, only bootstrap intervals are available.
#'
#' Appropriate sample quantities are used as starting values for all marginal parameters, regardless of the level of measurement. Details are provided in the package vignette.
#'
#' @param data a matrix of scores. Each row corresponds to a unit, each column a coder. The columns must be named appropriately so that the correct copula correlation matrix can be constructed. See \code{\link{build.R}} for details regarding column naming.
#' @param level the level of measurement, one of \code{"nominal"}, \code{"ordinal"}, \code{"interval"}, or \code{"ratio"}.
#' @param confint the method for computing confidence intervals, one of \code{"none"}, \code{"bootstrap"}, or \code{"asymptotic"}.
#' @param verbose logical; if TRUE, various messages are printed to the console.
#' @param control a list of control parameters.
#'    \describe{
#'        \item{\code{bootit}}{the size of the (parametric) bootstrap sample. This applies when \code{confint = "bootstrap"}, or when \code{confint = "asymptotic"} and \code{level = "nominal"} or \code{level = "ordinal"}. Defaults to 1,000.}
#'        \item{dist}{when \code{level = "interval"}, one of \code{"gaussian"}, \code{"laplace"}, \code{"t"}, \code{"gamma"}, or \code{"empirical"}.}
#'        \item{nodes}{the desired number of nodes in the cluster.}
#'        \item{parallel}{logical; if TRUE (the default), bootstrapping is done in parallel.}
#'        \item{type}{one of the supported cluster types for \code{\link[parallel]{makeCluster}}. Defaults to \code{"SOCK"}.}
#' }
#'
#' @return Function \code{sklars.omega} returns an object of class \code{"sklarsomega"}, which is a list comprising the following elements.
#'         \item{AIC}{the value of AIC for the fit, if \code{level = "interval"} and \code{dist != "empirical"}.}
#'         \item{BIC}{the value of BIC for the fit, if \code{level = "interval"} and \code{dist != "empirical"}.}
#'         \item{boot.sample}{when applicable, the bootstrap sample.}
#'         \item{call}{the matched call.}
#'         \item{coefficients}{a named vector of parameter estimates.}
#'         \item{confint}{the value of argument \code{confint}.}
#'         \item{control}{the list of control parameters.}
#'         \item{convergence}{unless optimization failed, the value of \code{convergence} returned by \code{\link{optim}}.}
#'         \item{cov.hat}{if \code{confint = "asymptotic"}, the estimate of the covariance matrix of the parameter estimator.}
#'         \item{data}{the matrix of scores, perhaps altered to remove rows (units) containing fewer than two scores.}
#'         \item{iter}{if optimization converged, the value of \code{iter} returned by \code{\link{optim}}.}
#'         \item{level}{the level of measurement.}
#'         \item{message}{the value of \code{massage} returned by \code{\link{optim}}.}
#'         \item{method}{the approach to inference, one of \code{"CML"}, \code{"DT"}, \code{"ML"}, or \code{"SMP"} (semiparametric).}
#'         \item{mpar}{the number of marginal parameters.}
#'         \item{npar}{the total number of parameters.}
#'         \item{R}{the initial value of the copula correlation matrix.}
#'         \item{R.hat}{the estimated value of the copula correlation matrix.}
#'         \item{residuals}{the residuals.}
#'         \item{root.R.hat}{a square root of the estimated copula correlation matrix. This is used for simulation and to compute the residuals.}
#'         \item{value}{the minimum of the log objective function.}
#'         \item{verbose}{the value of argument \code{verbose}.}
#'         \item{y}{the scores as a vector, perhaps altered to remove rows (units) containing fewer than two scores.}
#'
#' @references
#' Hughes, J. (2018). Sklar's Omega: A Gaussian copula-based framework for assessing agreement. \emph{ArXiv e-prints}, March.
#' @references
#' Krippendorff, K. (2013). Computing Krippendorff's alpha-reliability. Technical report, University of Pennsylvania.
#' @references
#' Nissi, M. J., Mortazavi, S., Hughes, J., Morgan, P., and Ellermann, J. (2015). T2* relaxation time of acetabular and femoral cartilage with and without intra-articular Gd-DTPA2 in patients with femoroacetabular impingement. \emph{American Journal of Roentgenology}, \bold{204}(6), W695.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # The following data were presented in Krippendorff (2013).
#'
#' data = matrix(c(1,2,3,3,2,1,4,1,2,NA,NA,NA,
#'                 1,2,3,3,2,2,4,1,2,5,NA,3,
#'                 NA,3,3,3,2,3,4,2,2,5,1,NA,
#'                 1,2,3,3,2,4,4,1,2,5,1,NA), 12, 4)
#' colnames(data) = c("c.1.1", "c.2.1", "c.3.1", "c.4.1")
#' data
#'
#' # Compute asymptotic confidence intervals. Since the distributional transform
#' # approximation is used, the asymptotic covariance matrix has a sandwich form.
#' # We use a bootstrap estimator of the "filling" (the variance of the score).
#' # Obtain a bootstrap sample of size 1,000. Do the bootstrap in parallel, using
#' # all but one of the available cores. Since we set 'verbose' equal to TRUE, a
#' # progress bar is displayed during the bootstrap.
#'
#' fit = sklars.omega(data, level = "nominal", confint = "asymptotic", verbose = TRUE,
#'                    control = list(bootit = 1000, parallel = TRUE,
#'                    nodes = parallel::detectCores() - 1))
#' summary(fit)
#' }
#' # Fit a subset of the cartilage data, assuming a Laplace marginal distribution. Compute
#' # confidence intervals in the usual ML way (observed information matrix).
#'
#' data(cartilage)
#' data = as.matrix(cartilage)[1:100, ]
#' colnames(data) = c("c.1.1", "c.2.1")
#' fit = sklars.omega(data, level = "interval", confint = "asymptotic",
#'                    control = list(dist = "laplace"))
#' summary(fit)
#' vcov(fit)
#'
#' # Now assume a t marginal distribution.
#'
#' fit = sklars.omega(data, level = "interval", confint = "asymptotic",
#'                    control = list(dist = "t"))
#' summary(fit)
#' vcov(fit)

sklars.omega = function(data, level = c("nominal", "ordinal", "interval", "ratio"), confint = c("none", "bootstrap", "asymptotic"),
                        verbose = FALSE, control = list())
{
    call = match.call()
    cnames = check.colnames(data)
    if (missing(data) || ! is.matrix(data) || ! is.numeric(data))
        stop("You must supply a numeric data matrix.")
    if (! cnames$success)
        stop("Your data matrix must have appropriately named columns. Check column(s) ", paste(cnames$cols, collapse = ", "), ".")
    level = match.arg(level)
    confint = match.arg(confint)
    if (! is.logical(verbose) || length(verbose) > 1)
        stop("'verbose' must be a logical value.")
    if (! is.list(control))
        stop("'control' must be a list.")
    if (level == "interval" && confint == "asymptotic" && ! is.null(control$dist) && control$dist == "empirical")
    {
        if (verbose)
        {
            cat("\nParameter 'confint' may not be equal to \"asymptotic\" if control parameter 'dist' is equal to \"empirical\".")
            cat(" Setting 'confint' equal to \"bootstrap\".\n")
        }
        confint = "bootstrap"
    }
    control = sklarsomega.control(level, confint, verbose, control)
    temp = build.R(data)
    R = temp$R
    onames = temp$onames
    rm(temp)
    block.sizes = apply(data, 1, function(row) { sum(! is.na(row)) })
    data = data[block.sizes > 1, ]
    y = as.vector(t(data))
    y = y[! is.na(y)]
    if (level %in% c("nominal", "ordinal"))
    {
        if (length(unique(y)) < 5)
            method = "CML"
        else
            method = "DT"
    }
    else if (level == "interval")
    {
        if (control$dist != "empirical")
            method = "ML"
        else
        {
            method = "SMP"
            F.hat = ecdf(y)
        }
    }
    else if (level == "ratio")
        method = "ML"
    result = list()
    class(result) = "sklarsomega"
    m = length(unique(R[R != 0 & R != 1]))
    hessian = ifelse(confint != "none", TRUE, FALSE)
    if (method == "CML")
    {
        C = length(unique(y))
        p.init = table(y)
        p.init = p.init / sum(p.init)
        N = neighbor.list(R)
        K = matrix(0, 4, 2)
        fit = try(optim(par = c(rep(0.5, m), p.init), fn = objective.CML, y = y, R = R, N = N, K = K,
                        method = "L-BFGS-B", lower = c(rep(0.001, m), rep(0.001, C)), upper = c(rep(0.999, m), rep(0.999, C))), silent = TRUE)
        if (! inherits(fit, "try-error"))
        {
            p.hat = fit$par[-c(1:m)]
            p.hat = p.hat / sum(p.hat)
            fit$par[-c(1:m)] = p.hat
            if (hessian)
                fit$hessian = optimHess(fit$par, objective.CML, y = y, R = R, N = N, K = K)
            cnames = paste0("p", 1:C)
        }
    }
    else if (method == "DT")
    {
        C = length(unique(y))
        p.init = table(y)
        p.init = p.init / sum(p.init)
        fit = try(optim(par = c(rep(0.5, m), p.init), fn = objective.DT, y = y, R = R, method = "L-BFGS-B",
                        lower = c(rep(0.001, m), rep(0.001, C)), upper = c(rep(0.999, m), rep(0.999, C))), silent = TRUE)
        if (! inherits(fit, "try-error"))
        {
            p.hat = fit$par[-c(1:m)]
            p.hat = p.hat / sum(p.hat)
            fit$par[-c(1:m)] = p.hat
            if (hessian)
                fit$hessian = optimHess(fit$par, objective.DT, y = y, R = R)
            cnames = paste0("p", 1:C)
        }
    }
    else if (method == "ML")
    {
        if (control$dist == "beta")
        {
            temp1 = mean(y)
            temp2 = var(y)
            alpha0 = temp1 * ((temp1 * (1 - temp1)) / temp2 - 1)
            beta0 = (1 - temp1) * ((temp1 * (1 - temp1)) / temp2 - 1)
            init = c(alpha0, beta0)
            lower = c(rep(0.001, m), rep(0.001, 2))
            upper = c(rep(0.999, m), rep(Inf, 2))
            cnames = c("alpha", "beta")
        }
        else if (control$dist == "t")
        {
            init = c(median(abs(y - median(y))), median(y))
            lower = c(rep(0.001, m), 0.001, -Inf)
            upper = c(rep(0.999, m), rep(Inf, 2))
            cnames = c("nu", "mu")
        }
        else if (control$dist == "gamma")
        {
            temp1 = mean(y)
            temp2 = var(y)
            alpha0 = temp1^2 / temp2
            beta0 = temp1 / temp2
            init = c(alpha0, beta0)
            lower = c(rep(0.001, m), rep(0.001, 2))
            upper = c(rep(0.999, m), rep(Inf, 2))
            cnames = c("alpha", "beta")
        }
        else
        {
            init = c(mean(y), sd(y))
            lower = c(rep(0.001, m), -Inf, 0.001)
            upper = c(rep(0.999, m), rep(Inf, 2))
            cnames = c("mu", "sigma")
        }
        fit = try(optim(par = c(rep(0.5, m), init), fn = objective.ML, y = y, R = R, dist = control$dist, hessian = hessian,
                        method = "L-BFGS-B", lower = lower, upper = upper), silent = TRUE)
    }
    else
    {
        fit = try(optim(par = rep(0.5, m), fn = objective.ML, y = y, R = R, dist = control$dist, F.hat = F.hat, hessian = FALSE,
                        method = "L-BFGS-B", lower = rep(0.001, m), upper = rep(0.999, m)), silent = TRUE)
        cnames = NULL
    }
    if (inherits(fit, "try-error"))
    {
        warning("Optimization failed.")
        result$message = fit[1]
        result$call = call
        class(result) = c(class(result), "try-error")
        return(result)
    }
    if (fit$convergence != 0)
    {
        warning("Optimization failed to converge.")
        result$convergence = fit$convergence
        result$message = fit$message
        result$call = call
        return(result)
    }
    npar = length(fit$par)
    result$level = level
    result$verbose = verbose
    result$npar = npar
    result$coefficients = fit$par
    names(result$coefficients) = c(onames, cnames)
    result$mpar = length(cnames)
    result$confint = confint
    result$convergence = fit$convergence
    result$message = fit$message
    result$R = R
    vals = sort(unique(R[R != 0 & R != 1]))
    result$R.hat = R
    for (j in 1:m)
        result$R.hat[result$R == vals[j]] = result$coefficients[j]
    temp = svd(result$R.hat)
    result$root.R.hat = temp$u %*% diag(sqrt(temp$d), nrow(result$R))
    rm(temp)
    result$value = fit$value
    result$iter = fit$counts[1]
    if (level == "interval" && method != "SMP")
    {
        result$AIC = 2 * result$value + 2 * npar
        result$BIC = 2 * result$value + log(length(y)) * npar
    }
    result$data = data
    result$y = y
    rest = result$coefficients[(m + 1):npar]
    u = switch(control$dist,
               gaussian = pnorm(y, mean = rest[1], sd = rest[2]),
               laplace = plaplace(y, location = rest[1], scale = rest[2]),
               t = pt(y, df = rest[1], ncp = rest[2]),
               beta = pbeta(y, shape1 = rest[1], shape2 = rest[2]),
               gamma = pgamma(y, shape = rest[1], rate = rest[2]),
               empirical = F.hat(y),
               categorical = (sklarsomega::pcat(y, rest) + sklarsomega::pcat(y - 1, rest)) / 2)
    z = qnorm(u)
    z = ifelse(z == -Inf, qnorm(0.0001), z)
    z = ifelse(z == Inf, qnorm(0.9999), z)
    residuals = try(solve(result$root.R.hat) %*% z, silent = TRUE)
    if (! inherits(residuals, "try-error"))
        result$residuals = as.vector(residuals)
    result$method = method
    result$call = call
    result$control = control
    if (confint != "none")
    {
        if (confint == "bootstrap")
            result$boot.sample = bootstrap.helper(result)
        else
        {
            bread = try(solve(fit$hessian), silent = TRUE)
            if (inherits(bread, "try-error"))
            {
                warning("Inversion of the Hessian matrix failed.")
                result$cov.hat = bread[1]
            }
            else if (method == "ML")
                result$cov.hat = bread
            else
            {
                meat = sandwich.helper(result)
                result$cov.hat = bread %*% meat %*% bread
            }
        }
    }
    result
}


#' Print a summary of a Sklar's Omega fit.
#'
#' @details Unless optimization of the objective function failed, this function prints a summary of the fit. First, the value of the objective function at its maximum is displayed, along with the number of iterations required to find the maximum. Then the values of the control parameters (defaults and/or values supplied in the call) are printed. Then a table of estimates is shown. If applicable, the table includes confidence intervals. Finally, the values of \code{\link{AIC}} and BIC are displayed (if the scores were interval and inference was parametric).
#'
#' @param object an object of class \code{sklarsomega}, the result of a call to \code{\link{sklars.omega}}.
#' @param alpha the significance level for the confidence intervals. The default is 0.05.
#' @param digits the number of significant digits to display. The default is 4.
#' @param \dots additional arguments.
#'
#' @seealso \code{\link{sklars.omega}}
#'
#' @method summary sklarsomega
#'
#' @references
#' Nissi, M. J., Mortazavi, S., Hughes, J., Morgan, P., and Ellermann, J. (2015). T2* relaxation time of acetabular and femoral cartilage with and without intra-articular Gd-DTPA2 in patients with femoroacetabular impingement. \emph{American Journal of Roentgenology}, \bold{204}(6), W695.
#'
#' @export
#'
#' @examples
#' # Fit a subset of the cartilage data, assuming a Laplace marginal distribution. Compute
#' # confidence intervals in the usual ML way (observed information matrix).
#'
#' data(cartilage)
#' data = as.matrix(cartilage)[1:100, ]
#' colnames(data) = c("c.1.1", "c.2.1")
#' fit = sklars.omega(data, level = "interval", confint = "asymptotic",
#'                    control = list(dist = "laplace"))
#' summary(fit)
#' vcov(fit)

summary.sklarsomega = function(object, alpha = 0.05, digits = 4, ...)
{
    cat("\nCall:\n\n")
    print(object$call)
    cat("\nConvergence:\n")
    if (is.null(object$convergence))
    {
        cat("\nOptimization failed.\n")
        return(invisible())
    }
    else if (object$convergence != 0)
    {
        cat("\nOptimization failed =>", object$message, "\n")
        return(invisible())
    }
    else
        cat("\nOptimization converged at", signif(-object$value, digits = digits), "after", object$iter, "iterations.\n")
    cat("\nControl parameters:\n")
    if (length(object$control) > 0)
    {
        control.table = cbind(unlist(c(object$control, "")))
        colnames(control.table) = ""
        print(control.table, quote = FALSE)
    }
    else
        print("\nNone specified.\n")
    npar = object$npar
    if (object$confint == "none")
        confint = matrix(rep(NA, 2 * npar), ncol = 2)
    else
    {
        boot.sample = object$boot.sample
        if (! is.null(boot.sample))
        {
            se = apply(boot.sample, 2, sd)
            scale = qnorm(1 - alpha / 2)
            coef = object$coef
            confint = cbind(coef - scale * se, coef + scale * se)
        }
        else
        {
            cov.hat = object$cov.hat
            if (! is.null(cov.hat))
            {
                se = sqrt(diag(cov.hat))
                scale = qnorm(1 - alpha / 2)
                coef = object$coef
                confint = cbind(coef - scale * se, coef + scale * se)
            }
            else
                confint = matrix(rep(NA, 2 * npar), ncol = 2)
        }
    }
    coef.table = cbind(object$coef, confint)
    colnames(coef.table) = c("Estimate", "Lower", "Upper")
    rownames(coef.table) = names(object$coef)
    cat("Coefficients:\n\n")
    print(signif(coef.table, digits = digits))
    if (! is.null(object$AIC))
        cat("\nAIC:", signif(object$AIC, digits), "\nBIC:", signif(object$BIC, digits), "\n")
    cat("\n")
}


#' Extract model residuals.
#'
#' @details Although our simulation studies suggest that residuals are not terribly useful in this context, we provide residuals nonetheless. Said residuals are computed by first applying the probability integral transform, then applying the inverse probability integral transform, then pre-multiplying by the inverse of the square root of the (fitted) copula correlation matrix. For nominal or ordinal scores, the distributional transform approximation is used.
#'
#' @param object an object of class \code{sklarsomega}, typically the result of a call to \code{\link{sklars.omega}}.
#' @param \dots additional arguments.
#'
#' @return A vector of residuals.
#'
#' @seealso \code{\link{sklars.omega}}
#'
#' @method residuals sklarsomega
#'
#' @references
#' Nissi, M. J., Mortazavi, S., Hughes, J., Morgan, P., and Ellermann, J. (2015). T2* relaxation time of acetabular and femoral cartilage with and without intra-articular Gd-DTPA2 in patients with femoroacetabular impingement. \emph{American Journal of Roentgenology}, \bold{204}(6), W695.
#'
#' @export
#'
#' @examples
#' # Fit a subset of the cartilage data, assuming a Laplace marginal distribution.
#'
#' data(cartilage)
#' data = as.matrix(cartilage)[1:100, ]
#' colnames(data) = c("c.1.1", "c.2.1")
#' fit = sklars.omega(data, level = "interval", control = list(dist = "laplace"))
#' summary(fit)
#' res = residuals(fit)
#' qqnorm(res, pch = 20)
#' abline(0, 1, col = "red", lwd = 2)

residuals.sklarsomega = function(object, ...)
{
    if (! is.null(object$residuals))
        return(object$residuals)
    else
        stop("Fit object does not contain a vector of residuals.")
}


#' Compute an estimated covariance matrix for a Sklar's Omega fit.
#'
#' @details See the package vignette for detailed information regarding covariance estimation for Sklar's Omega.
#'
#' @param object a fitted model object.
#' @param ... additional arguments.
#
#' @return A matrix of estimated variances and covariances for the parameter estimator. This should have row and column names corresponding to the parameter names given by the \code{\link{coef}} method. Note that a call to this function will result in an error if \code{\link{sklars.omega}} was called with argument \code{confint} equal to \code{"none"}, or if optimization failed.
#'
#' @method vcov sklarsomega
#'
#' @references
#' Nissi, M. J., Mortazavi, S., Hughes, J., Morgan, P., and Ellermann, J. (2015). T2* relaxation time of acetabular and femoral cartilage with and without intra-articular Gd-DTPA2 in patients with femoroacetabular impingement. \emph{American Journal of Roentgenology}, \bold{204}(6), W695.
#'
#' @export
#'
#' @examples
#' # Fit a subset of the cartilage data, assuming a Laplace marginal distribution. Compute
#' # confidence intervals in the usual ML way (observed information matrix).
#'
#' data(cartilage)
#' data = as.matrix(cartilage)[1:100, ]
#' colnames(data) = c("c.1.1", "c.2.1")
#' fit = sklars.omega(data, level = "interval", confint = "asymptotic",
#'                    control = list(dist = "laplace"))
#' summary(fit)
#' vcov(fit)

vcov.sklarsomega = function(object, ...)
{
    if (is.null(object$confint))
        stop("Fit object does not contain field 'confint'\n.")
    if (object$confint == "none")
        stop("Function sklars.omega was called with 'confint' equal to \"none\".")
    if (! is.null(object$cov.hat))
        cov.hat = object$cov.hat
    else
        cov.hat = cov(object$boot.sample, use = "complete.obs")
    rownames(cov.hat) = colnames(cov.hat) = names(object$coef)
    cov.hat
}


#' Simulate a Sklar's Omega dataset(s).
#'
#' @details This function simulates one or more responses distributed according to the fitted model.
#'
#' @param object a fitted model object.
#' @param nsim number of datasets to simulate. Defaults to 1.
#' @param seed either \code{NULL} or an integer that will be used in a call to \code{\link{set.seed}} before simulating the response vector(s). If set, the value is saved as the \code{"seed"} attribute of the returned value. The default (\code{NULL}) will not change the random generator state, and \code{\link{.Random.seed}} will be returned as the \code{"seed"} attribute.
#' @param ... additional arguments.
#
#' @return A data frame having \code{nsim} columns, each of which contains a simulated response vector. Said data frame has a \code{"seed"} attribute, which takes the value of the \code{seed} argument or the value of \code{\link{.Random.seed}}.
#'
#' @method simulate sklarsomega
#'
#' @export
#'
#' @examples
#' # The following data were presented in Krippendorff (2013).
#'
#' data = matrix(c(1,2,3,3,2,1,4,1,2,NA,NA,NA,
#'                 1,2,3,3,2,2,4,1,2,5,NA,3,
#'                 NA,3,3,3,2,3,4,2,2,5,1,NA,
#'                 1,2,3,3,2,4,4,1,2,5,1,NA), 12, 4)
#' colnames(data) = c("c.1.1", "c.2.1", "c.3.1", "c.4.1")
#' fit = sklars.omega(data, level = "nominal", confint = "none")
#' summary(fit)
#'
#' # Simulate three datasets from the fitted model, and then
#' # display the first dataset in matrix form.
#'
#' sim = simulate(fit, nsim = 3, seed = 12)
#' data.sim = t(fit$data)
#' data.sim[! is.na(data.sim)] = sim[, 1]
#' data.sim = t(data.sim)
#' data.sim

simulate.sklarsomega = function(object, nsim = 1, seed = NULL, ...)
{
    n = nrow(object$R)
    psi = object$coef[(object$npar - object$mpar + 1):object$npar]
    data = NULL
    if (! is.null(seed))
       set.seed(seed)
    else
       seed = .Random.seed
    j = 1
    while (j <= nsim)
    {
        z = as.numeric(object$root.R.hat %*% rnorm(n))
        u = pnorm(z)
        y = switch(object$control$dist,
                   gaussian = qnorm(u, mean = psi[1], sd = psi[2]),
                   laplace = qlaplace(u, location = psi[1], scale = psi[2]),
                   t = qt(u, df = psi[1], ncp = psi[2]),
                   beta = qbeta(u, shape1 = psi[1], shape2 = psi[2]),
                   gamma = qgamma(u, shape = psi[1], rate = psi[2]),
                   empirical = quantile(object$y, u, type = 8),
                   categorical = sklarsomega::qcat(u, p = psi))
        if (object$control$dist != "categorical" ||
            (object$control$dist == "categorical" && length(unique(y)) == object$mpar))
        {
            data = cbind(data, as.numeric(y))
            j = j + 1
        }
    }
    data = as.data.frame(data)
    colnames(data) = paste0("sim_", 1:ncol(data))
    attr(data, "seed") = seed
    data
}


#' Return the number of observations for a Sklar's Omega fit.
#'
#' @details This function extracts the number of observations from a model fit, and is principally intended to be used in computing information criteria.
#'
#' @param object a fitted model object.
#' @param ... additional arguments.
#
#' @return An integer.
#'
#' @seealso \code{\link{AIC}}
#'
#' @method nobs sklarsomega
#'
#' @export

nobs.sklarsomega = function(object, ...)
{
    length(object$y)
}


#' Return the maximum of the Sklar's Omega log objective function.
#'
#' @details This function extracts the maximum value of the log objective function from a model fit, and is principally intended to be used in computing information criteria.
#'
#' @param object a fitted model object.
#' @param ... additional arguments.
#
#' @return This function returns an object of class \code{logLik}. This is a number with at least one attribute, \code{"df"} (degrees of freedom), giving the number of estimated parameters in the model.
#'
#' @seealso \code{\link{AIC}}
#'
#' @method logLik sklarsomega
#'
#' @export

logLik.sklarsomega = function(object, ...)
{
    result = -object$value
    attr(result, "df") = object$npar
    class(result) = "logLik"
    result
}


#' Compute DFBETAs for units and/or coders.
#'
#' @details This function computes DFBETAS for one or more units and/or one or more coders.
#'
#' @param model a fitted model object.
#' @param units a vector of integers. A DFBETA will be computed for each of the corresponding units.
#' @param coders a vector of integers. A DFBETA will be computed for each of the corresponding coders.
#' @param ... additional arguments.
#
#' @return A list comprising at most two elements.
#'         \item{dfbeta.units}{a matrix, the columns of which contain DFBETAS for the units specified via argument \code{units}.}
#'         \item{dfbeta.coders}{a matrix, the columns of which contain DFBETAS for the coders specified via argument \code{coders}.}
#'
#' @method influence sklarsomega
#'
#' @references
#' Young, D. S. (2017). \emph{Handbook of Regression Methods}. CRC Press.
#' @references
#' Krippendorff, K. (2013). Computing Krippendorff's alpha-reliability. Technical report, University of Pennsylvania.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # The following data were presented in Krippendorff (2013).
#'
#' data = matrix(c(1,2,3,3,2,1,4,1,2,NA,NA,NA,
#'                 1,2,3,3,2,2,4,1,2,5,NA,3,
#'                 NA,3,3,3,2,3,4,2,2,5,1,NA,
#'                 1,2,3,3,2,4,4,1,2,5,1,NA), 12, 4)
#' colnames(data) = c("c.1.1", "c.2.1", "c.3.1", "c.4.1")
#' fit = sklars.omega(data, level = "nominal", confint = "none")
#' summary(fit)
#' (inf = influence(fit, units = c(6, 11), coders = c(2, 3)))
#' }

influence.sklarsomega = function(model, units, coders, ...)
{
    dfbeta.units = NULL
    if (! missing(units))
    {
        dfbeta.units = matrix(NA, length(units), model$npar)
        k = 1
        for (j in units)
        {
            fit = try(suppressWarnings(sklars.omega(model$data[-j, ], level = model$level, confint = "none", control = model$control)), silent = TRUE)
            if (inherits(fit, "try-error") || fit$convergence != 0)
                warning(fit$message)
            else
                dfbeta.units[k, ] = model$coef - fit$coef
            k = k + 1
        }
        colnames(dfbeta.units) = names(model$coef)
        rownames(dfbeta.units) = units
    }
    dfbeta.coders = NULL
    if (! missing(coders))
    {
        dfbeta.coders = matrix(NA, length(coders), model$npar)
        k = 1
        for (j in coders)
        {
            fit = try(suppressWarnings(sklars.omega(model$data[, -j], level = model$level, confint = "none", control = model$control)), silent = TRUE)
            if (inherits(fit, "try-error") || fit$convergence != 0)
                warning(fit$message)
            else
                dfbeta.coders[k, ] = model$coef - fit$coef
            k = k + 1
        }
        colnames(dfbeta.coders) = names(model$coef)
        rownames(dfbeta.coders) = coders
    }
    result = list()
    if (! is.null(dfbeta.units))
        result$dfbeta.units = dfbeta.units
    if (! is.null(dfbeta.coders))
        result$dfbeta.coders = dfbeta.coders
    result        
}


#' Produce a Bland-Altman plot.
#'
#' @details This function produces rather customizable Bland-Altman plots, using the \code{\link{plot}} and \code{\link{abline}} functions. The former is used to create the scatter plot. The latter is used to display the confidence band.
#'
#' @param x the first vector of outcomes.
#' @param y the second vector of outcomes.
#' @param pch the plotting character. Defaults to 20, a bullet.
#' @param col the foreground color. Defaults to black.
#' @param bg the background color. Defaults to black.
#' @param main a title for the plot. Defaults to no title, i.e., \code{""}.
#' @param xlab a title for the x axis. Defaults to \code{"Mean"}.
#' @param ylab a title for the y axis. Defaults to \code{"Difference"}.
#' @param lwd1 the line width for the scatter plot. Defaults to 1.
#' @param lwd2 the line width for the confidence band. Defaults to 1.
#' @param cex scaling factor for the scatter plot. Defaults to 1.
#' @param lcol line color for the confidence band. Defaults to black.
#
#' @references
#' Altman, D. G. and Bland, J. M. (1983). Measurement in medicine: The analysis of method comparison studies. \emph{The Statistician}, 307--317.
#' @references
#' Nissi, M. J., Mortazavi, S., Hughes, J., Morgan, P., and Ellermann, J. (2015). T2* relaxation time of acetabular and femoral cartilage with and without intra-articular Gd-DTPA2 in patients with femoroacetabular impingement. \emph{American Journal of Roentgenology}, \bold{204}(6), W695.
#'
#' @export
#'
#' @examples
#' # Reproduce the plot from Figure 4 of the package vignette.
#'
#' data(cartilage)
#' baplot(cartilage$pre, cartilage$post, pch = 21, col = "navy", bg = "darkorange", lwd1 = 2,
#'        lwd2 = 2, lcol = "navy")

baplot = function(x, y, pch = 20, col = "black", bg = "black", main = "", xlab = "Mean",
                  ylab = "Difference", lwd1 = 1, lwd2 = 1, cex = 1, lcol = "black")
{
    bamean = (x + y) / 2
    badiff = y - x
    plot(badiff ~ bamean, pch = pch, xlab = xlab, ylab = ylab, lwd = lwd1, main = main, col = col, bg = bg, cex = cex)
    abline(h = c(mean(badiff), mean(badiff) + 2 * sd(badiff), mean(badiff) - 2 * sd(badiff)), lty = 2, lwd = lwd2, col = lcol)
} 


