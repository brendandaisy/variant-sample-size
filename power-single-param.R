library(tidyverse)
library(cowplot)

source("src/growth-model.R")
source("src/nopt-single-param.R")

## Code to test whether the framework matches the pairwise normals implementation
## in the two region case:

# r <- c(0.12, 0.07)
# q0 <- 1e-4
# tmax <- 60
# n_prop <- c(1/2, 1/2)
# 
# nopt_multiple_regions(tmax, q0, r, 2, H0_true=FALSE, power=0.8)
# 
# vc <- diag(map_dbl(r, ~var_single_param(tmax, q0, .)))
# sd <- sqrt(sum(vc/(30*n_prop)))
# crit <- qnorm(0.025, 0, sd)
# pnorm(crit, r[2]-r[1], sd) # yay!

# nopt_power_right <- function(tmax, q0, r, power=0.8, level=0.05) {
#     s2 <- var_single_param(tmax, q0, r)
#     Z1 <- qnorm(1-level) # notice one sided this time
#     Z2 <- qnorm(power)
#     s2 / r^2 * (Z1 + Z2)^2
# }
# TODO: use me!
wald_stat <- function(tmax, q0, r, M) {
    vc <- diag(map2_dbl(q0, r, ~var_single_param(tmax, .x, .y)))
    C <- cbind(rep(-1, M-1), diag(1, M-1))
    W <- solve(C %*% vc %*% t(C))
    t(C %*% r) %*% W %*% (C %*% r) # test statistic with n=1 samples
}

nopt_multiple_regions <- function(tmax, q0, r_alt, M, n_prop=rep(1/M, M), power=0.8, level=0.05) {
    if (!isTRUE(all.equal(sum(n_prop), 1)))
        stop("n_prop must be a valid probability vector")
    if (length(r_alt) == 1) {
        r_alt <- rep(r_alt, M)
    }
    if (length(q0) == 1)
        q0 <- rep(q0, M)
    if (length(r_alt) != M | length(q0) != M | length(n_prop) != M)
        stop("incompatible number of regions across r, q0, n_prop")

    vc <- diag(map2_dbl(q0, r_alt, ~var_single_param(tmax, .x, .y)))
    C <- cbind(rep(-1, M-1), diag(1, M-1))
    W <- solve(C %*% (vc/n_prop) %*% t(C))
    tstat <- t(C %*% r_alt) %*% W %*% (C %*% r_alt) # test statistic with n=1 samples
    crit_t <- qchisq(1-level, M-1) # critical value we need tstat to exceed to reject null
    
    n <- 1
    # to update the tstat, multiply by n cause can pull out the inverse C * Sigma C'
    target_prob <- pchisq(crit_t, M-1, ncp=n*tstat, lower.tail=FALSE)
    while (target_prob < power) {
        n <- 2*n
        target_prob <- pchisq(crit_t, M-1, ncp=n*tstat, lower.tail=FALSE)
    }
    
    n_rem <- n - n/2
    n <- n - floor(n_rem/2) # half point between n (too big) and n/2 (too small)
    while (n_rem > 1) {
        target_prob <- pchisq(crit_t, M-1, ncp=n*tstat, lower.tail=FALSE)
        if (target_prob > power) # n too big
            n <- n - n_rem
        else # n too small
            n <- n + n_rem
        n_rem <- floor(n_rem/2)
    }
    
    return(n)
    # if (H0_true)
    #     ret <- nopt_multiple_regions_null(tmax, q0, r, M, n_prop, power, level)
    # else
    #     ret <- nopt_multiple_regions_alt(tmax, q0, r, M, n_prop, power, level)
    # 
    # return(ret)
}

# `level` is our desired type I error!
power_multiple_regions <- function(tmax, q0, r_alt, M, n, level=0.05) {
    if (length(r_alt) == 1) {
        r_alt <- rep(r_alt, M)
    }
    if (length(q0) == 1)
        q0 <- rep(q0, M)
    if (length(r_alt) != M | length(q0) != M | length(n) != M)
        stop("incompatible number of regions across r_alt, q0, n_prop")
    
    vcn <- diag(map2_dbl(q0, r_alt, ~var_single_param(tmax, .x, .y))) / n
    # ^ notice already adjusting for sampling size / proportions
    C <- cbind(rep(-1, M-1), diag(1, M-1))
    W <- solve(C %*% vcn %*% t(C))
    tstat <- t(C %*% r_alt) %*% W %*% (C %*% r_alt) # test statistic with reported n
    
    crit_t <- qchisq(1-level, M-1) # critical value we need tstat to exceed to reject null
    beta <- pchisq(crit_t, M-1, ncp=tstat)
    1 - beta
}

# nopt_multiple_regions_alt <- function(tmax, q0, r, M, n_prop, power, level) {
#     vc <- diag(map2_dbl(q0, r, ~var_single_param(tmax, .x, .y)))
#     C <- cbind(rep(-1, M-1), diag(1, M-1))
#     W <- solve(C %*% (vc/n_prop) %*% t(C))
#     tstat <- t(C %*% r) %*% W %*% (C %*% r) # test statistic with n=1 samples
#     crit_t <- qchisq(1-level, M-1) # critical value we need tstat to exceed to reject null
#     
#     n <- 1
#     # to update the tstat, multiply by n cause can pull out the inverse C * Sigma C'
#     target_prob <- pchisq(crit_t, M-1, ncp=n*tstat, lower.tail=FALSE)
#     while (target_prob < power) {
#         n <- 2*n
#         target_prob <- pchisq(crit_t, M-1, ncp=n*tstat, lower.tail=FALSE)
#     }
#     
#     n_rem <- n - n/2
#     n <- n - floor(n_rem/2) # half point between n (too big) and n/2 (too small)
#     while (n_rem > 1) {
#         target_prob <- pchisq(crit_t, M-1, ncp=n*tstat, lower.tail=FALSE)
#         if (target_prob > power) # n too big
#             n <- n - n_rem
#         else # n too small
#             n <- n + n_rem
#         n_rem <- floor(n_rem/2)
#     }
#     
#     return(n)
# }

###Two? example scenarios which we will present for the MA data:

## Determine whether a variant, newly established in region A, has an inflated
# growth advantage compared to other regions
M <- 3 # num counties
n_prop <- rep(1, M) / M # assume even sampling

## Assume an Omicron-like r is detected in A, and what to know if it ends up 
## being as bad in other regions
## Specify the assumed r in A, which has reached 5%, and is has now reached 0.1% in others
r0 <- 0.12
r <- c(r0, rep(r0-r0*0.15, M-1))
q0 <- c(0.05, rep(1e-3, M-1))
# q0 <- c(0.05, 1e-3, 0.01, 1e-4)
tmax <- 30

nopt_multiple_regions(tmax, q0, r, M, n_prop) / M # samples per region

###
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

#  Now for power with multiple regions--------------------------------------------

r <- c(0.12, 0.07)
q0 <- 1e-4
tmax <- 60
n_prop <- c(1/2, 1/2)

nopt_multiple_regions(tmax, q0, r, 2, H0_true=FALSE, power=0.8)

vc <- diag(map_dbl(r, ~var_single_param(tmax, q0, .)))
sd <- sqrt(sum(vc/(30*n_prop)))
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

