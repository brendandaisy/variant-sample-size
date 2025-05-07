library(tidyverse)
library(cowplot)

source("src/growth-model.R")
source("src/nopt-single-param.R")

# i.e. "right sided" hypothesis r > 0
# TODO: notice level is reversed from orig. too
nopt_power_right <- function(tmax, q0, r, power=0.8, level=0.05) {
    s2 <- var_single_param(tmax, q0, r)
    Z1 <- qnorm(1-level) # notice one sided this time
    Z2 <- qnorm(power)
    s2 / r^2 * (Z1 + Z2)^2
}

q0 <- 1e-4
r <- 0.11
tmax <- 40

nopt_power_right(tmax, q0, r, power=0.8)
nopt_single_param(tmax, q0, r, error=0.05)

nopt_multiple_regions <- function(tmax, q0, r, M, n_prop=rep(1/M, M), 
                                  H0_true=TRUE, power=0.8, level=0.05) {
    if (!isTRUE(all.equal(sum(n_prop), 1)))
        stop("n_prop must be a valid probability vector")
    if (length(r) == 1) {
        r <- rep(r, M)
    }
    if (length(q0) == 1)
        q0 <- rep(q0, M)
    if (H0_true)
        ret <- nopt_multiple_regions_null(tmax, q0, r, M, n_prop, power, level)
    else
        ret <- nopt_multiple_regions_alt(tmax, q0, r, M, n_prop, power, level)
    
    return(ret)
}

nopt_multiple_regions_alt <- function(tmax, q0, r, M, n_prop, power, level) {
    vc <- diag(map2_dbl(q0, r, ~var_single_param(tmax, .x, .y)))
    C <- cbind(rep(-1, M-1), diag(1, M-1))
    W <- solve(C %*% (vc/n_prop) %*% t(C))
    tstat <- t(C %*% r) %*% W %*% (C %*% r) # test statistic with n=1 samples
    crit_t <- qchisq(1-level, M-1) # critical value we need tstat to exceed to reject null
    
    n <- 1
    target_prob <- pchisq(crit_t, M-1, ncp=n*tstat, lower.tail=FALSE)
    while (target_prob < power) {
        n <- 2*n
        target_prob <- pchisq(crit_t, M-1, ncp=n*tstat, lower.tail=FALSE)
    }
    
    n_rem <- n - n/2
    n <- n - floor(n_rem/2) # half point between n (too big) and n/2 (too small)
    while (n_rem > 1) {
        target_prob <- pchisq(crit_t, M-1, ncp=n*tstat, lower.tail=FALSE)
        n_rem <- floor(n_rem/2)
        if (target_prob > power) # n too big
            n <- n - n_rem
        else # n too small
            n <- n + n_rem
    }
    return(n)
}

r <- c(0.11, 0.11)
q0 <- 1e-4
tmax <- 80
n <- 1
n_prop <- c(0.01, 0.99)

nopt_2reg_res <- tibble(np1=seq(0.01, 0.99, 0.01), np2=1-np1) |> 
    rowwise() |> 
    mutate(nopt=nopt_multiple_regions(tmax, q0, r, 2, c(np1, np2))) |> 
    ungroup()

nopt_2reg_min <- slice_min(nopt_2reg_res, nopt, with_ties=FALSE)

ggplot(nopt_2reg_res, aes(np1, nopt)) +
    geom_line() +
    geom_point(data=nopt_2reg_min, col=pal_discrete[6], size=2) +
    scale_x_continuous(sec.axis=sec_axis(~1-., name="proportion region 2")) +
    labs(x="proportion region 1", y="number of sequences per day")
    # coord_cartesian(ylim=c(0, 2000))

nopt_multiple_regions(tmax, q0, r, 2, n_prop)

vc <- diag(map_dbl(r, ~var_single_param(tmax, q0, .)))
C <- matrix(c(-1, 1), nrow=1)

wald_info <- solve(C %*% (vc/n_prop) %*% t(C))

tstat <- 5*t(C %*% r) %*% wald_info %*% (C %*% r)


pf(tstat, 1, 259-2, lower.tail=FALSE)
pchisq(tstat, 1, lower.tail=FALSE) 

sd <- sqrt(sum(vc/(n_fact*n_rate)))
2*pnorm(r[2]-r[1], 0, sd) # matches for M=2 version!

#  Now for power with multiple regions--------------------------------------------

r <- c(0.12, 0.07, 0.07)
q0 <- 1e-4
tmax <- 60
n_prop <- c(1/4, 1/2, 1/4)

nopt_multiple_regions(tmax, q0, r, 3, n_prop, H0_true=FALSE, power=0.8)

sd <- sqrt(sum(vc/(305*n_prop)))
crit <- qnorm(0.025, 0, sd)
pnorm(crit, r[2]-r[1], sd) # yay!

f_nopt <- function(q0, r, M, n_prop) nopt_multiple_regions(tmax, q0, r, M, n_prop, H0_true=FALSE, power=0.8)

nopt_res_mult_reg <- tibble(M=2:10) |> 
    rowwise() |> 
    mutate(
        r=list(c(0.12, rep(0.07, M-1))),
        q0=list(c(1e-4, rep(0.8, M-1))),
        n_prop_under=list(c(0.05, rep(0.95/(M-1), M-1))),
        n_prop_even=list(rep(1/M, M)),
        n_prop_over=list(c(0.85, rep(0.15/(M-1), M-1))),
        nopt_under=f_nopt(q0, r, M, n_prop_under),
        nopt_even=f_nopt(q0, r, M, n_prop_even),
        nopt_over=f_nopt(q0, r, M, n_prop_over),
    ) |> 
    ungroup()

nopt_res_mult_reg |> 
    pivot_longer(contains("nopt")) |> 
    ggplot(aes(M, value*tmax/M, col=name)) +
    geom_line(linewidth=1.05) +
    geom_point(size=2.1) +
    # scale_y_continuous(breaks=0:3*1000, limits=c(0, NA)) +
    scale_color_manual(
        values=c(pal_discrete[c(2, 3, 5)]), 
        labels=c("even distribution", "region A oversampled", "region A undersampled")
    ) +
    labs(x="number of regions", y="total sequences per region", col=NULL) +
    theme_half_open() +
    theme(legend.position=c(0.46, 0.8), plot.margin=unit(c(1, 16, 1, 1), "mm"))

ggsave("figs/power-mult-reg-over-under.pdf", width=4.5, height=3.7)

## compare power to simulated expectation TODO: full version
library(MASS, exclude=c("select"))

rsim <- mvrnorm(100000, r, vc/(n_fact*n_rate))
colnames(rsim) <- c("r1", "r2")

rsim_rej <- as_tibble(rsim) |> 
    mutate(diff=r2-r1, n=n()) |> 
    filter(diff < crit)

nrow(rsim_rej) / 100000

# 
# ggplot(as_tibble(rsim), aes(r1, r2)) +
#     stat_ellipse() +
#     geom_abline(intercept=0, slope=1)

