library(tidyverse)
# library(matrixStats) # for logSumExp

# part one: derive _simulation_ based probability of success for a logistic prevalence curve for their
# prescribed n. If they do not match we know something is up
growth_model <- function(ts, q0, r) {
    a <- 1/q0 - 1
    1 / (1 + a*exp(-r * ts))
}

gm_proflik_optim <- function(y, n, q0=1e-4) {
    if (all(y == 0))
        return(tibble(t0=length(y)-1, r=0)) # basically a definition
    rmax <- optimize(
        function(r) gm_log_lik(r, y, n, q0),
        interval=c(0, 3),
        maximum=TRUE
    )$maximum
    t0max <- gm_t0_est(rmax, y, n, q0)
    return(tibble(t0=t0max, r=rmax))
}

# assumes length of y and ts are the same
gm_log_lik <- function(r, y, n, q0=1e-4) {
    ts <- seq_along(y)-1
    Tmax <- last(ts)
    t0 <- gm_t0_est(r, y, n, q0)
    q <- c(rep(0, t0), growth_model(0:(Tmax-t0-1), q0, r))
    sum(dbinom(y, n, q, log=TRUE))
    # if (last(ts) < t0) {
    #     if (any(y > 0))
    #         return(-Inf)
    #     return(0)
    # }
    # if (r <= 0)
    #     return(-Inf)
    # ts_growth <- ts[ts >= t0] - t0
    # q <- c(rep(0, sum(ts < t0)), growth_model(ts_growth, q0, r))
    # sum(dbinom(y, n, q, log=TRUE))
}

gm_t0_est <- function(r, y, n, q0=1e-4) {
    ts <- seq_along(y)-1
    Tmax <- last(ts)
    q <- growth_model(ts, q0, r)
    cond_lik <- map_dbl(ts, \(t0) {
        qshift <- c(rep(0, t0), q[1:(Tmax-t0+1)])
        # print(paste0(t0, " --- ", qshift))
        sum(dbinom(y, n, qshift, log=TRUE))
    })
    ts[which.max(cond_lik)]
}

gm_marg_log_lik <- function(r, y, n, q0=1e-4) {
    ts <- seq_along(y)-1
    q <- growth_model(ts, q0, r)
    
    ts_growth <- ts[ts >= t0] - t0
    sum(dbinom(y, n, q, log=TRUE))
}

simulate_data <- function(tmax, n, t0_true, q0_true, r_true, nsim=10000) {
    ts <- 0:tmax # timesteps we assume sequencing will take place
    ts_growth_true <- 0:(tmax-t0_true-1)
    q_true <- c(rep(0, t0_true), growth_model(ts_growth_true, q0_true, r_true))
    
    y <- map(1:nsim, \(rep) rbinom(tmax+1, n, q_true))
    return(list(y=y, n=n, t0=t0_true, q0=q0_true, r=r_true, tmax=tmax, qmax=last(q_true)))
}

t <- seq(0, 50, 1)
qt <- growth_model(t, 1/10000, 0.2)
y <- rbinom(0:50, 50, qt)

ggplot(tibble(t, qt, y), aes(t, qt)) +
    geom_line() +
    geom_point(aes(y=y/50), col="purple") +
    labs(x="Days of sequencing", y="Latent/observed prevalence") +
    theme_bw() +
    # theme(plot.background=element_blank(), panel.border=element_blank(), panel.background=element_blank())
    theme(rect=element_rect(fill="transparent"))

voc_detect_prob <- function(tmax, q, n, nsim=10000) {
    sims <- map_lgl(1:nsim, ~any(rbinom(tmax, n, q[1:tmax]) > 0))
    mean(sims)
}

voc_growth_ci_mc <- function(sim_true, nsim=length(ytrue), level=0.95) {
    ysub <- sample(sim_true$y, nsim) # shuffle/downsample true data
    mle <- map_dfr(ysub, \(y) gm_proflik_optim(y, sim_true$n, sim_true$q0))
    
    quants <- c((1 - level)/2, (1 + level)/2)
    ci <- quantile(mle$r, quants)
    return(list(
        y=ysub,
        mle=mle,
        ci=ci,
        error=max(abs(ci - sim_true$r)),
        cv=sd(mle$r) / sim_true$r,
        cv_obs=sd(mle$r) / mean(mle$r)
    ))
}

dom_time <- function(q0, r) {
    a <- 1/q0 - 1
    log(a) / r
}

voc_dom_time <- function(tmax, q0_true, r_true, n, nsim=10000, level=0.95) {
    ts <- 0:tmax # timesteps we assume sequencing will take place
    q_true <- growth_model(ts, q0_true, r_true)
    dom_true <- dom_time(q0_true, r_true)
    
    dom_times <- map_dbl(1:nsim, \(rep) {
        y <- rbinom(tmax+1, n, q_true)
        # ML values for q0 and r
        mlp <- optim(
            c(q0_true, r_true), function(p) gm_log_lik(p, y, n, ts), 
            method="L-BFGS-B",
            lower=c(1e-7, 1e-6), upper=c(0.2, 0.8), # TODO: fix problem of high r = error (q=1 but y=0)
            control=list(fnscale=-1)
        )$par
        
        dom_time(mlp[1], mlp[2])
    })
    
    quants <- c((1 - level)/2, (1 + level)/2)
    ci <- quantile(dom_times, quants)
    return(list(
        mle=dom_times,
        ci=ci,
        error=max(abs(ci - dom_true))
    ))
}

find_sample_size <- function(criterion, q, tmax, nsim=10000) {
    # todo not really useful since ur gonna plot as func of sample size anyway
}

# TODO plot as function of n and tmax - how many samples needed for various windows/prevalences?
voc_detect_prob(30, qt, 200)

sim_true <- simulate_data(45, 100, 5, 1e-4, 0.1, nsim=10000)

res_growth <- voc_growth_ci_mc(sim_true, nsim=2000)

ggplot(res_growth$mle, aes(r, t0)) +
    geom_density_2d_filled() +
    annotate("point", x=0.1, y=5, col="tomato2", shape=8)

ggplot(res_growth$mle, aes(r)) +
    geom_histogram() +
    geom_vline(xintercept=0.1, col="tomato2", linetype="dashed")

res_growth$error

###
design_grid <- expand_grid(tmax=c(30, 60, 5), n=c(100, 250, 400))

design_grid |> 
    mutate(err=map2_dbl())
map2_dfr(design_grid$tmax, design_grid$n, \(tmax, n) {
    sim_true <- simulate_data(tmax, n, 5, 1e-4, 0.1, nsim=10000)
    res <- voc_growth_ci_mc(sim_true, nsim=2000)
})

# how many samples are needed, *as a function of the final prevalence/days of observation*,
# to have a 95% chance of an error less than 0.05 for the growth rate, assuming a true growth rate 
# of 0.1?

# (this is definitely where an analytical formula would be useful)

function(tmax, error=0.05, level=0.95) { # tmax should be more than t0=5
    nprev <- 0
    ncurr <- 100
    nmax <- 10000
    sat_prev <- FALSE
    while ((ncurr < nmax) & (ncurr - nprev) > 10) { # TODO check this
        n <- ncurr
        sim_true <- simulate_data(tmax, n, 5, 1e-4, 0.1, nsim=5000)
        res <- voc_growth_ci_mc(sim_true, nsim=500, level=level)
        
        if (res$error <= error) {
            ncurr <- ceiling((ncurr-nprev)/2)
        }
        else {
            ncurr <- 2*ncurr
        }
        nprev <- n
    }
    
}

###

res_dom <- voc_dom_time(60, 1e-4, 0.1, 50, nsim=5000)

ggplot(tibble(x=res_dom$mle), aes(x)) +
    geom_density(adjust=1, bounds=c(0, 250)) +
    geom_vline(xintercept=dom_time(1e-4, 0.1), linewidth=1.04, linetype="dashed", col="purple") +
    # coord_cartesian(xlim=c(0, 200)) +
    labs(x="estimated dom. time", y="density") +
    theme_bw()
