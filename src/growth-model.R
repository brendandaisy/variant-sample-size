library(tidyverse)

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
    if (t0 < 0) { # TODO this case relatively untested
        ts <- seq(-t0, -t0+tmax, 1)
        return(1 / (1 + a*exp(-r * ts)))
    }
    ts <- seq((1 - t0)%%1, tmax-t0, 1) # assumes Tmax is an integer
    c(rep(0, ceiling(t0)), 1 / (1 + a*exp(-r * ts)))
}

library(microbenchmark)
microbenchmark(
    gm_old=growth_model(30, 5, 1e-4, 0.1),
    gm_new=growth_model2(30, 5, 1e-4, 0.1)
)

final_prevalence <- function(tmax, t0, q0, r) {
    a <- 1/q0 - 1
    1 / (1 + a*exp(-r * (tmax-t0)))
}

simulate_data <- function(tmax, n, t0_true, q0_true, r_true, nsim=10000) {
    q_true <- growth_model2(tmax, t0_true, q0_true, r_true)
    
    y <- map(1:nsim, \(rep) rbinom(length(q_true), n, q_true))
    return(list(y=y, n=n, t0=t0_true, q0=q0_true, r=r_true, tmax=tmax, qmax=last(q_true)))
}

# Maximum likelihood estimation---------------------------------------------------
gm2_optim <- function(y, n, tmax, q0=1e-4, r_ubound=3) {
    if (all(y == 0))
        return(tibble(t0=length(y)-1, r=0)) # basically a definition
    mle <- optim(
        c(runif(1, 1, tmax/4), 0.1 + rnorm(1, 0, 0.05)),
        \(theta) gm2_log_lik(theta, y, n, tmax, q0),
        method="BFGS",
        hessian=FALSE, 
        control=list(fnscale=-1, ndeps=c(1e-4, 1e-4))
        # lower=c(0, 0), upper=c(tmax-1, r_ubound)
    )
    return(tibble(t0=mle$par[1], r=mle$par[2]))
}

gm2_log_lik <- function(theta, y, n, tmax, q0) {
    # if (theta[1] < 0 | theta[1] >= tmax)
    #     return(-Inf)
    q <- growth_model2(tmax, theta[1], q0, theta[2])
    # print(paste0(theta[2], " -- ", last(q)))
    # not super safe but a convenient way to defend against q=0 or 1
    # max(sum(dbinom(y, n, q, log=TRUE)), -1e12)
    sum(dbinom(y, n, q, log=TRUE))
}

# sim <- simulate_data(60, 100, 4.9, 1e-4, 0.12, nsim=100)
# gm2_optim(sim$y[[4]], 100, 60)

## functions for discrete t0

gm_log_lik <- function(r, y, n, q0=1e-4) {
    tmax <- length(y) - 1
    t0 <- gm_t0_est(r, y, n, q0)
    q <- growth_model(tmax, t0, q0, r)
    sum(dbinom(y, n, q, log=TRUE))
}

gm_t0_est <- function(r, y, n, q0=1e-4) {
    tmax <- length(y) - 1
    t0s <- 0:tmax
    cond_lik <- map_dbl(t0s, \(t0) {
        q <- growth_model(tmax, t0, q0, r)
        # print(paste0(t0, " --- ", qshift))
        sum(dbinom(y, n, q, log=TRUE))
    })
    t0s[which.max(cond_lik)]
}

gm_proflik_optim <- function(y, n, q0=1e-4, r_int=c(0, 3)) {
    if (all(y == 0))
        return(tibble(t0=length(y)-1, r=0)) # basically a definition
    rmax <- optimize(
        function(r) gm_log_lik(r, y, n, q0),
        interval=r_int,
        maximum=TRUE
    )$maximum
    t0max <- gm_t0_est(rmax, y, n, q0)
    return(tibble(t0=t0max, r=rmax))
}