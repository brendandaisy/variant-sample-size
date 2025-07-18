library(tidyverse)
library(cowplot)

source("src/growth-model.R")
source("src/nopt-single-param.R")

# Comparing asymptotic and Monte Carlo error as n increases-----------------------
mc_err_single_param <- function(sim, r_bound=c(-2, 2), level=0.95) {
    r_mle <- with(sim, {
        map_dbl(y, \(.y) {
            # if (all(.y == 0)) { # incorrect
            #     return(0)
            # }
            fun <- function(r) sum(dbinom(.y, n, growth_model2(tmax, t0, q0, r), log=TRUE))
            optimize(fun, r_bound, maximum=TRUE)$maximum
        })
    })
    
    quants <- c((1 - level)/2, (1 + level)/2)
    ci <- quantile(r_mle, quants)
    max(abs(ci - sim$r))
}

q0 <- 1e-4
# r <- c(0.07, 0.12) * 7
r <- 0.12*7
tmax <- 8

# approx_err <- map2_dfr(params$qmax, params$r, \(qT, r) {
#     tmax <- find_tmax(qT, 0, q0, r)
#     Z <- -qnorm((1 - 0.95)/2)
#     ns <- seq(10, 100, 1)
#     map_dfr(ns, \(n) {
#         tibble(
#             tmax=tmax, qmax=qT, r=r, n=n,
#             err=Z * sqrt(var_single_param(n, tmax, 0, q0, r))
#         )
#     })
# })
# 
# ggplot(approx_err, aes(n, err, col=factor(r), linetype=factor(qmax))) +
#     geom_line() +
#     scale_linewidth_manual(values=c("1", "31")) +
#     # scale_color_manual(values=c("darkorange1", "#35A0ABFF")) +
#     labs(
#         x="sequences per day",
#         y="error (95% CI)"
#         # col="selection\ncoefficient",
#     ) +
#     theme_half_open()

Z <- -qnorm((1 - 0.95)/2)
err_form <- expand_grid(r, ns=seq(1, 100, 1)) |> 
    rowwise() |> 
    mutate(
        n=ns,
        err_raw={
            vars <- var_single_param(tmax, q0, r)
            Z * sqrt(vars/n)
        },
        err=err_raw / r * 100
    )

err_mc <- expand_grid(tmax, ns=c(3, seq(5, 90, 5))) |> 
    rowwise() |> 
    mutate(
        n=ns,
        sim=list(simulate_data(tmax, n, 0, q0, r, nsim=5000)),
        err_mc_raw=mc_err_single_param(sim, r_bound=c(-0.25, 0.25)),
        err_mc=err_mc_raw / r * 100
    )

# for presentation:
sd1 <- sqrt(var_single_param(20, tmax[2], 0, q0, r))

sim <- simulate_data(tmax[2], 20, 0, q0, r, nsim=5000)

sim_bad <- err_mc |>
    filter(tmax == 76, n == 3) |> 
    pull(sim)

r_mle <- with(sim_bad[[1]], {
    map_dbl(y, \(.y) {
        fun <- function(r) sum(dbinom(.y, 20, growth_model2(tmax, t0, q0, r), log=TRUE))
        optimize(fun, c(-0.25, 0.25), maximum=TRUE)$maximum
    })
})

ggplot() +
    geom_histogram(aes(r_mle, after_stat(density)), tibble(r_mle), bins=50, col="gray20", fill="gray70") +
    # stat_function(fun=~dnorm(., r, sd1), n=500, linewidth=1.2, col="#35A0ABFF") +
    # xlim(0.04, 0.1) +
    labs(x="estimated growth rate", y=NULL) +
    theme_half_open()

# ggsave("figs/example-asymp-approx.pdf", width=3.6, height=3.6)

# panel B:
ggplot(err_form, aes(n, err, col=r)) +
    # geom_hline(yintercept=25, col="gray70", linewidth=1.2, linetype="dashed") +
    geom_hline(yintercept=15, col="gray70", linewidth=1.2, linetype="dashed") +
    # geom_point(aes(y=err_mc), err_mc) +
    # geom_line(aes(y=err_mc), err_mc, linetype="dashed") +
    geom_line(linewidth=1.2, col="#F6B245") +
    scale_y_continuous(labels=scales::label_percent(scale=1)) +
    # scale_color_manual(values=c("darkorange1", "#35A0ABFF")) +
    labs(x="sequences per week", y="error", col="final\nprevalence") +
    # coord_cartesian(ylim=c(0, 100)) +
    theme_half_open()

ggsave("figs/error-vs-n.pdf", width=3.8, height=2.5)

ggsave("figs/error-asymp-vs-mc-2.pdf", width=5, height=3.5)

#  Prevalence uncertainty and sample size for increasing t------------------------
r <- 0.12 # daily sequencing for this example (cleaner curve)
q0 <- 5e-3
tmax <- find_tmax(0.99, q0, r)
n <- 6

p1 <- ggplot(var_prev_t) +
    stat_function(
        fun=\(x) dnorm(x, r, sqrt(var_single_param(30, q0, r)/n)), 
        col=pal_discrete[1], linewidth=1.2
        # geom="area", fill=pal_discrete[1]
    ) +
    stat_function(
        fun=\(x) dnorm(x, r, sqrt(var_single_param(55, q0, r)/n)), 
        col=pal_discrete[6], linewidth=1.2
        # geom="area", fill=pal_discrete[1]
    ) +
    xlim(r+c(-0.04, 0.04)) +
    labs(x="growth advantage", y="density") +
    theme_half_open()

var_prev_t <- tibble(t=1:tmax, q=growth_model(tmax, q0, r)[-1]) |>
    rowwise() |>
    mutate(
        vr=var_single_param(t, q0, r)/n,
        vprev=vr * d_prev_single_param(t, q0, r)^2,
        ymin=qnorm(0.025, q, sqrt(vprev)), 
        ymax=qnorm(0.975, q, sqrt(vprev)),
        dprev=d_prev_single_param(t, q0, r)
    ) |> 
    ungroup()

mdprev <- max(var_prev_t$dprev)

p2 <- ggplot(var_prev_t, aes(t)) +
    geom_vline(xintercept=30, linetype="dashed", linewidth=1.2, col=pal_discrete[1]) +
    geom_vline(xintercept=55, linetype="dashed", linewidth=1.2, col=pal_discrete[6]) +
    geom_line(aes(y=dprev/mdprev), col="gray50", linewidth=1.1, linetype="dotted") +
    geom_ribbon(aes(ymin=ymin, ymax=ymax), fill=pal_discrete[3], col=NA, alpha=0.65) +
    geom_line(aes(y=q), col=pal_discrete[3], linewidth=1.2) +
    scale_y_continuous(sec.axis=sec_axis(transform=~.*mdprev, name="derivative")) +
    # scale_x_continuous(breaks=c(0, 5, 10)) +
    labs(x="days", y="prevalence") +
    coord_cartesian(ylim=c(0, 1)) +
    theme_half_open()

plot_grid(p1, p2, nrow=1, rel_widths=c(0.4, 0.6), labels=c("A", "B"))

ggsave("figs/delta-method-ex.pdf", width=8, height=3.4)
