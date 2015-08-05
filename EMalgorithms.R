# TODO: implement the sweep operator to deal with NAs.
ECME <- function (T.sample, tol = 1e-5, nu.t = 100, maxiter = 500) {
# Author: G.filion
# Date: June 24, 2009
# Modified EM algorithm for T distributed variables. Returns the
# estimated mean, variance and degree of freedom of the distribution.

    E.step <- function (mu.t, var.t, nu.t) {

        w <<- ( nu.t + 1 ) / ( nu.t + (T.sample - mu.t)^2 / var.t )
        S.tau <<- sum(w)
        S.tau.X <<- sum(w * T.sample)
        S.tau.XX <<- sum(w * T.sample^2)

    }

    n <- length(T.sample)

    mu.t <- mean(T.sample)
    var.t <- var(T.sample)

    mu.t.plus.1 = mu.t + 1
    var.t.plus.1 = var.t + 1

    for (i in 1:maxiter) {
        E.step (mu.t, var.t, nu.t)

        # CM-step 1

        mu.t.plus.1 = S.tau.X / S.tau
        var.t.plus.1 = ( S.tau.XX - S.tau.X^2 / S.tau ) / n

        E.step (mu.t.plus.1, var.t.plus.1, nu.t)

        # CM-step 2

        # Estimate nu.
        RHS <- sum(log(w) - w) / n + 1
        no.upper.bound <- TRUE
        lower.bound <- 0
        # Find an upper bound.
        nu <- nu.t
        while (no.upper.bound) {
            if (digamma(nu/2) - log(nu/2) - digamma((nu + 1)/2) + log((nu + 1)/2) < RHS) {
                lower.bound <- nu
                nu <- 2 * nu
            }
            else {
                no.upper.bound <- FALSE
                upper.bound <- nu
            }
        }
        # The degree of freedom nu is estimated with a precision of 0.01
        while (upper.bound - lower.bound > 0.01) {
            nu <- (upper.bound + lower.bound) / 2
            if (digamma(nu/2) - log(nu/2) - digamma((nu + 1)/2) + log((nu + 1)/2) < RHS)
                lower.bound <- nu
            else
                upper.bound <- nu
        }
        nu.t <- round(nu, 1)
        cat(paste(mu.t, "\n"))

        if (abs(mu.t.plus.1 - mu.t) && abs(var.t.plus.1 - var.t) < tol)
            break

        mu.t <- mu.t.plus.1
        var.t <- var.t.plus.1

    }

    return (list(mu = mu.t.plus.1, var = var.t.plus.1, nu = nu.t))

}


EM.negative.binomial <- function (x, zero.censored = FALSE, tol = 1e-5) {

# Author: G. Filion
# Date: July 24, 2009
# EM as described in An EM algorithm for estimating negative binomial parameters
# K. Adamidis, Austral & New Zealand J. Statist. 41(2), 1999, 213-221
# Parameters:
#   x: either the sample, or a table of the observations.
#   zero.censored: T if zero is not observed.
#   tol: stop when max{|theta.t-theta.t+1|, |lamda.t-lambda.t+1|} < tol.

    if (!is.table(x))
        x <- table(x)

    obs <- as.integer(names(x))
    counts <- as.integer(x)
    n <- sum(counts)

    counts <- counts[obs != 0]
    obs <- obs[obs != 0]
        
    # Define a constant used in the E-step
    list.1 <- lapply(as.list(obs), function(x) seq(from=0, to=x-1))


    # Find the initial estimates by moment estimation
    x.bar   <- sum(counts*obs) / n
    sigma.2 <- sum(counts*obs^2) / n - x.bar^2
    new.theta   <- 1 - x.bar / sigma.2
    alpha   <- -1/log(1-new.theta)          #alias
    new.lambda  <- x.bar*(1-new.theta)/new.theta / alpha
    
    E.step <- function () {
        
        alpha <- - 1/log(1-theta)       # alias
        E.Z   <<- -((1-theta)/theta - alpha)
        E.Mi  <<- alpha*lambda*
            sapply(lapply(list.1, function(x) 1/(x + alpha*lambda)), sum)

    }

    lambda <- new.lambda + tol + 1
    theta  <- new.theta + tol + 1

    i = 0
    while (any(abs(c(theta, lambda) - c(new.theta, new.lambda)) > tol)) {

        i = i+1
        cat(paste("Cycle", i, "\n"))

        theta <- new.theta
        lambda <- new.lambda

        E.step()

        # M.step
        Cst = ifelse(zero.censored, 1-exp(-lambda), 1)
        new.theta  <- sum(counts*(obs-E.Mi)) / sum(counts*(obs+(E.Z-1)*E.Mi))
        new.lambda <- sum(counts*E.Mi) / n * Cst

    }

    return(c(alpha = -new.lambda / log(1-new.theta), theta = new.theta))

}
