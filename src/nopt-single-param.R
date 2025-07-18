# putting functions for error, power, and all single param methods in one place

var_single_param <- function(tmax, q0, r) {
    a <- 1/q0 - 1
    ts <- 0:tmax
    Ir <- sum(ts^2 * (a*exp(-r*ts)) / (1 + a*exp(-r*ts))^2)
    return(1 / Ir)
}

nopt_single_param <- function(tmax, q0, r, error=0.05, level=0.05) {
    vr <- var_single_param(tmax, q0, r)
    Z <- qnorm(1 - level/2)
    return((Z/error)^2 * vr)
}

nopt_g_single_param <- function(t, q0, r, dgdr, error=0.05, level=0.05, tmax=t) { # generally would assume we only observe up to t as well
    vr <- var_single_param(tmax, q0, r)
    vg <- vr * dgdr(t, q0, r)^2
    Z <- qnorm(1 - level/2)
    return((Z/error)^2 * vg)
}

d_prev_single_param <- function(t, q0, r) {
    a <- 1/q0 - 1
    eterm <- exp(-r * t)
    denom <- (1 + a*eterm)^2
    dqdr <- a * t * eterm / denom
    return(dqdr)
}

d_dom_single_param <- function(t, q0, r) { # TODO unused t somewhat awkward
    a <- 1/q0 - 1
    -log(a) / r^2
}

d_tmax_single_param <- function(q, t, q0, r) {
    a <- 1/q0 - 1
    -log(-(q - 1)/(a*q))/r^2
}

nopt_prev_classic <- function(t, q0, r, error=0.05, level=0.05) {
    Z <- qnorm(1 - level/2)
    qtrue <- final_prevalence(t, q0, r)
    (Z/error)^2 * qtrue * (1 - qtrue)
}

