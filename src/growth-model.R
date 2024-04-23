library(tidyverse)

# TODO: maybe t0 can be continuous?

# 3-4 times as fast as the ifelse version:
growth_model <- function(tmax, t0, q0, r) {
    a <- 1/q0 - 1
    c(rep(0, t0), 1 / (1 + a*exp(-r * 0:(tmax-t0))))
}

# a continuous version cause why not
# currently is like 20 times slower
growth_model2 <- function(tmax, t0, q0, r) {
    a <- 1/q0 - 1
    # if (isTRUE(all.equal(t0, floor(t0))))
    #     return(c(rep(0, t0), 1 / (1 + a*exp(-r * 0:(tmax-t0)))))
    ts <- seq((1 - t0)%%1, tmax-t0, 1) # assumes Tmax is an integer
    c(rep(0, ceiling(t0)), 1 / (1 + a*exp(-r * ts)))
}

# microbenchmark(
#     gm_old=growth_model(30, 5, 1e-4, 0.1),
#     gm_new=growth_model2(30, 5, 1e-4, 0.1)
# )

final_prevalence <- function(tmax, t0, q0, r) {
    last(growth_model(tmax, t0, q0, r))
}

simulate_data <- function(tmax, n, t0_true, q0_true, r_true, nsim=10000) {
    q_true <- growth_model(tmax, t0_true, q0_true, r_true)
    
    y <- map(1:nsim, \(rep) rbinom(tmax+1, n, q_true))
    return(list(y=y, n=n, t0=t0_true, q0=q0_true, r=r_true, tmax=tmax, qmax=last(q_true)))
}