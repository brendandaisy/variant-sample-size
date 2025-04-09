library(tidyverse)

growth_model <- function(tmax, q0, r) {
    a <- 1/q0 - 1
    1 / (1 + a*exp(-r * 0:tmax))
}

# library(microbenchmark)
# microbenchmark(
#     gm_old=growth_model(30, 5, 1e-4, 0.1),
#     gm_new=growth_model2(30, 5, 1e-4, 0.1)
# )

final_prevalence <- function(tmax, q0, r) {
    a <- 1/q0 - 1
    1 / (1 + a*exp(-r * tmax))
}

find_tmax <- function(qmax, q0, r) {
    a <- 1/q0 - 1
    ceiling(-log((1-qmax)/(a*qmax)) / r)
}

dom_time <- function(t0, q0, r) {
    # r < 0 and t0 != 0 untested
    stop("fix me please :(")
    a <- 1/q0 - 1
    stopifnot("r must be nonzero"=r != 0)
    ceiling((log(a) + r*t0) / r)
}

# find_min_r_tmax <- function(qmax, tmax, t0, q0) {
#     a <- 1/q0 - 1
# }

# simulate_data <- function(tmax, n, t0_true, q0_true, r_true, nsim=10000) {
#     q_true <- growth_model2(tmax, t0_true, q0_true, r_true)
#     
#     y <- map(1:nsim, \(rep) rbinom(length(q_true), n, q_true))
#     return(list(y=y, n=n, t0=t0_true, q0=q0_true, r=r_true, tmax=tmax, qmax=last(q_true)))
# }