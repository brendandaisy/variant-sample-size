# putting functions for error, power, and all single param methods in one place

var_single_param <- function(tmax, q0, r) {
    a <- 1/q0 - 1
    ts <- 0:tmax
    Ir <- sum(ts^2 * (a*exp(-r*ts)) / (1 + a*exp(-r*ts))^2)
    return(1 / Ir)
}

nopt_single_param <- function(tmax, q0, r, error=0.05, level=0.95) {
    vr <- var_single_param(tmax, q0, r)
    Z <- -qnorm((1 - level)/2)
    return((Z/error)^2 * vr)
}

