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

dom_time <- function(q0, r) {
    # r < 0 and t0 != 0 untested
    if (q0 == 0.5)
        return(0)
    if (r == 0)
        return(Inf)
    # if (q0 < 0.5 & r < 0)
    #     return(Inf)
    # if (q0 > 0.5 & r > 0)
    #     return(0)
    a <- 1/q0 - 1
    ceiling(log(a) / r)
}

info_per_t <- function(t, q0, r) {
    a <- 1/q0 - 1
    eterm <- exp(-r*t)
    (t^2*a*eterm) / (1 + a*eterm)^2
}

# solve for where the derivative of I_t(r) == 0 in t
# f is said derivative
info_root_t <- function(q0, r) {
    a <- 1/q0 - 1
    tmax <- find_tmax(0.99, q0, r)
    f <- function(t) (a*t*exp(r*t) * (a*(r*t + 2) + exp(r*t)*(2 - r*t)))/(a + exp(r*t))^3
    uniroot(f, c(1, tmax), tol=0.00001)
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