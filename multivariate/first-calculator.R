library(tidyverse)
library(furrr)
library(cowplot)

source("src/growth-model.R")

# TODO: may need to update for unknown t0?
# voc_detect_prob <- function(tmax, q, n, nsim=10000) {
#     sims <- map_lgl(1:nsim, ~any(rbinom(tmax, n, q[1:tmax]) > 0))
#     mean(sims)
# }

# `nsim` number of samples from sim to use. Default is to shuffle all of them. Set to NULL to not shuffle
# `...` control arguments passed to gm_proflik_optim
voc_growth_ci_mc <- function(sim, nsim=length(ytrue), level=0.95) {
    if (!is.null(nsim)) {
        y <- sample(sim$y, nsim) # shuffle/downsample true data
    } else
        y <- sim$y
    
    # mle <- furrr::future_map_dfr(y, \(y) gm_proflik_optim(y, sim$n, sim$q0, ...))
    mle <- map_dfr(
        y, 
        \(y) gm2_optim(y, sim$n, sim$tmax, sim$q0)
        # .options=furrr_options(seed=TRUE)
    )
    
    quants <- c((1 - level)/2, (1 + level)/2)
    ci <- quantile(mle$r, quants)
    return(list(
        y=y,
        mle=mle,
        ci=ci,
        error=max(abs(ci - sim$r)),
        cv=sd(mle$r) / sim$r,
        cv_obs=sd(mle$r) / mean(mle$r)
    ))
}

voc_fin_prev_ci_mc <- function(sim_true, nsim=length(ytrue), level=0.95, ...) {
    if (!is.null(nsim)) {
        y <- sample(sim_true$y, nsim) # shuffle/downsample true data
    } else
        y <- sim_true$y
    
    mle <- furrr::future_map_dfr(y, \(y) gm_proflik_optim(y, sim_true$n, sim_true$q0, ...))
    qmax_mle <- map2_dbl(mle$t0, mle$r, \(t0, r) final_prevalence(sim_true$tmax, t0, sim_true$q0, r))
    
    quants <- c((1 - level)/2, (1 + level)/2)
    ci <- quantile(qmax_mle, quants)
    return(list(
        y=y,
        mle=qmax_mle,
        ci=ci,
        error=max(abs(ci - sim_true$qmax)),
        cv=sd(qmax_mle) / sim_true$qmax,
        cv_obs=sd(qmax_mle) / mean(qmax_mle)
    ))
}

sim_true <- simulate_data(35, 1000, 5, 1e-4, 0.12, nsim=10000)

plan(multisession, workers=8)
res_growth <- voc_growth_ci_mc(sim_true, nsim=NULL, r_int=c(0, 1))

ggplot(res_growth$mle, aes(r, t0)) +
    geom_density_2d_filled() +
    annotate("point", x=0.1, y=5, col="tomato2", shape=8)

ggplot(res_growth$mle, aes(r)) +
    geom_histogram() +
    geom_vline(xintercept=0.1, col="tomato2", linetype="dashed")

res_growth$error
res_growth$cv

# r_true <- 0.1
# q0_true <- 1e-4
# t0_true <- 5
# 
# design_grid <- expand_grid(n=round(seq(50, 5000, length.out=15)), tmax=round(seq(30, 90, length.out=10)))
# 
# gr_err_grid <- map2_dfr(design_grid$n, design_grid$tmax, \(n, tmax) {
#     print(paste0("n=", n, ", tmax=", tmax))
#     # fin_prev <- last(growth_model(0:(tmax-t0_true-1), q0_true, r_true))
#     sim_true <- simulate_data(tmax, n, t0_true, q0_true, r_true, nsim=5000)
#     res_growth <- voc_growth_ci_mc(sim_true, nsim=500, r_int=c(0, 0.8))
#     tibble_row(
#         n, qmax=sim_true$qmax, tmax=tmax, 
#         err=res_growth$error, cv=res_growth$cv,
#         r_true=r_true, q0_true=q0_true, t0_true=t0_true
#     )
# })
# 
# gr_err_grid |> 
#     mutate(tmax=rep(seq(30, 90, length.out=10), times=15)) |> 
#     ggplot(aes(tmax, n, fill=err)) +
#     geom_raster() +
#     scale_fill_gradientn(
#         colors=c("cyan", "black", "red"), 
#         values=scales::rescale(c(0.03, 0.045, 0.05, 0.1, 0.4))
#     ) +
#     scale_x_continuous(expand=expansion()) +
#     scale_y_continuous(expand=expansion())

# optimal sample size stuff-------------------------------------------------------

# TODO 4/26: this needs to have an argument for the confidence level (not same as error!)
# `...` control arguments passed to gm_proflik_optim
find_sample_size <- function(sim_fun, obj_fun, error=0.05, nmax=10000, ntol=10, etol=1e-3) {
    # first check the case that n* >= nmax
    sim_true <- sim_fun(nmax)
    res <- obj_fun(sim_true)
    if (res$error > error)
        stop("nmax too small")
    
    # now we know n* < nmax and we proceed with a binary search
    rem <- floor(nmax/2)
    n <- nmax - rem
    while (rem > ntol) {
        sim_true <- sim_fun(n)
        res <- obj_fun(sim_true)
        rem <- floor(rem/2)
        print(paste0("n = ", n, ", err = ", res$error))
        if (abs(res$error - error) < etol) {
            return(list(nopt=n, res=res[-1]))
        }
        else if (res$error < error) { # n too large
            n <- n - rem
        }
        else { # n too small
            n <- n + rem
        }
    }
    # update the error for the final n value:
    sim_true <- sim_fun(n)
    res <- obj_fun(sim_true)
    
    return(list(nopt=n, res=res[-1]))
}
# find_sample_size <- function(sim_fun, error=0.05, level=0.95, nmax=10000, tol=10, ...) {
#     # first check the case that n* >= nmax
#     sim_true <- sim_fun(nmax)
#     res <- voc_growth_ci_mc(sim_true, nsim=NULL, level=level, ...)
#     if (res$error > error)
#         return("error - nmax too small")
#     
#     # now we know n* < nmax and we proceed with binary search
#     nleft <- 0
#     nright <- nmax
#     while ((nright-nleft) > tol) {
#         n <- floor((nright - nleft) / 2)
#         sim_true <- sim_fun(n)
#         res <- voc_growth_ci_mc(sim_true, nsim=NULL, level=level, ...)
#         
#         if (res$error < error) { # n too large
#             nright <- n - 1
#             print(nright)
#             print(nleft)
#             print(paste0("n = ", n, ", err = ", res$error))
#         }
#         else if (res$error > error) { # n too small
#             nleft <- n + 1
#             print(nright)
#             print(nleft)
#             print(paste0("n = ", n, ", err = ", res$error))
#         }
#         else
#             return(list(nopt=n, res=res[-1]))
#     }
#     return(list(nopt=nright, res=res[-1])) # err. on the side of n being a bit larger than n*
# }

r_alpha <- 0.07
r_delta <- 0.08
r_omicron <- 0.12

q0_true <- 1e-4
t0_true <- 5

prev_true <- imap_dfr(c(alpha=r_alpha, delta=r_delta, omicron=r_omicron), \(r, var) {
    tibble(variant=var, t=0:80, q=growth_model(80, t0_true, q0_true, r))
})

p1 <- ggplot(prev_true, aes(t, q, col=variant)) +
    geom_line() +
    scale_color_manual(values=c("tomato", "blue3", "green2"))

tmax <- seq(20, 80, 10)
nmax <- c(10000, 4000, 2000, 1000, 500, 200, 200)

# nopt_alpha <- numeric(length(tmax))
# nopt_delta <- numeric(length(tmax))
# nopt_omicron <- numeric(length(tmax))

obj_fun <- function(sim) voc_growth_ci_mc(sim, nsim=NULL, level=0.95, r_int=c(0, 1))
for (i in seq_along(tmax)) {
    sim_fun <- function(n) simulate_data(tmax[i], n, t0_true, q0_true, r_omicron, nsim=5000)
    nopt_res <- find_sample_size(sim_fun, obj_fun, nmax=nmax[i], ntol=5)
    nopt_omicron[i] <- nopt_res$nopt
}

nmax <- c(1000, 400, 200, 200, 200, 200, 400)

# nopt_alpha_qmax <- numeric(length(tmax))
# nopt_delta_qmax <- numeric(length(tmax))
# nopt_omicron_qmax <- numeric(length(tmax))

obj_fun <- function(sim) voc_fin_prev_ci_mc(sim, nsim=NULL, level=0.95, r_int=c(0, 1))
for (i in seq_along(tmax)) {
    sim_fun <- function(n) simulate_data(tmax[i], n, t0_true, q0_true, r_omicron, nsim=4000)
    nopt_res <- find_sample_size(sim_fun, obj_fun, nmax=nmax[i], error=0.01, ntol=5)
    nopt_omicron_qmax[i] <- nopt_res$nopt
}

nopt_res_alpha <- tibble(
    tmax, 
    nopt_growth=nopt_alpha, nopt_qmax=nopt_alpha_qmax, 
    qmax=map_dbl(tmax, ~final_prevalence(., t0_true, q0_true, r_alpha))
)

nopt_res_delta <- tibble(
    tmax, 
    nopt_growth=nopt_delta, nopt_qmax=nopt_delta_qmax, 
    qmax=map_dbl(tmax, ~final_prevalence(., t0_true, q0_true, r_delta))
)

nopt_res_omicron <- tibble(
    tmax, 
    nopt_growth=nopt_omicron, nopt_qmax=nopt_omicron_qmax, 
    qmax=map_dbl(tmax, ~final_prevalence(., t0_true, q0_true, r_omicron))
)

p4 <- nopt_res_omicron |> 
    pivot_longer(c(nopt_growth, nopt_qmax), values_to="nopt") |> 
    mutate(name=factor(name, labels=c("growth rate", "final prevalence"))) |> 
    ggplot(aes(tmax-5, nopt * tmax)) +
    geom_line(aes(linetype=name), col="green2") +
    geom_point(col="green2") +
    scale_x_continuous(
        breaks=tmax,
        sec.axis=dup_axis(name="final prevalence", labels=scales::scientific(nopt_res_omicron$qmax, 1))
    ) +
    labs(x="days since introduction", y="total sequences needed", linetype=NULL)

plot_grid(p1, p2, p3, p4, nrow=2, labels=c("A", "B", "C", "D"))
ggsave("figs/nopt-mc-qmax-err-pt01.pdf", width=9.5, height=5)


