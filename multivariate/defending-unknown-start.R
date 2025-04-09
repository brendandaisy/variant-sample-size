# analysis to justify our decision to consider the start time unknown
# a good supp fig 1
# 1. a figure illustrating the data generation process, motivating: "when did we start?"
# 2. we could analyse proportional data with 3 models: linear transform, simple logistic, logistic w/ unknown start

# supp fig 2: illustration of patching the logistic with unknown start, with q as func of t0

library(tidyverse)
library(cowplot)

source("src/growth-model.R")
# TODO: 11/17: for the presentation, most of the "full model" sample size calculations ended up at the top 
# of this script. Below that is a latest example of defending unknown start time that I am pretty happy
# with.

# 1/11: moving information matrix example plot here
r <- 0.12
q0 <- 1e-4
n <- 100
t0 <- 5.5
tmax <- find_tmax(0.1, t0, q0, r)

q_true <- growth_model2(tmax, t0, , r)

ex_dat <- tibble(
    q=q_true,
    y=rbinom(length(q), n, q),
    t=seq_along(y)-1
)

ll_grid <- expand_grid(t0=seq(-1.5, 6, 0.1), r=seq(0.2, 0.75, 0.02)) |> 
    mutate(loglik=map2_dbl(t0, r, ~gm2_log_lik(c(.x, .y), ex_dat$y, n, tmax, 1e-4)))

(im <- info_mat(tmax, t0, q0, r))

library(matrixcalc)
is.positive.definite(round(im, 2))

cov <- solve(round(im, 6))
mle_approx <- MASS::mvrnorm(10000, c(r, t0, 1/q0-1), cov)
colnames(mle_approx) <- c("r", "t0", "a")
mle_approx <- as_tibble(mle_approx)

Z * sqrt(cov[1, 1]/n)
Z * sqrt(solve(im[1:2, 1:2])[1, 1]/n) # thankfully doesn't impact inference for r very much
# err <- voc_growth_ci_im(n, 12, 2.5, 1e-4, r)

ggplot(mle_approx, aes(t0, a)) +
    # geom_raster(aes(fill=1.05^loglik), ll_grid, interpolate=TRUE) +
    stat_ellipse() +
    stat_ellipse(level=0.5, linetype="dotted") +
    # annotate("point", x=r, y=t0, size=1.8, shape=8) +
    # scale_fill_viridis_c(option="mako") +
    labs(x="selection coefficient", y="start time") +
    scale_x_continuous(expand=expansion()) +
    scale_y_continuous(expand=expansion()) +
    # coord_cartesian(ylim=c(NA, 6)) +
    theme_half_open() +
    theme(legend.position="none")

# 1/11: lesson learned (share with group!) t0 and q0 are highly collinear,
# but it doesn't appear to influence error for r that much so its probably ok

ggsave("figs/example-inf-mat.pdf", width=3, height=3)

### compare optimal sample size for full and simple models
r <- 0.12
q0 <- 1e-4
t0 <- 0

qmax <- seq(0.01, 0.75, length.out=60)

nopt_new_qmax <- tibble(qmax) |>
    rowwise() |>
    mutate(
        tmax=find_tmax(qmax, t0, q0, r),
        nopt_old=nopt_single_param(tmax, t0, q0, r, error=r*0.15),
        nopt_tot_old=nopt_old*tmax, 
        nopt=find_sample_size(tmax, t0, q0, r, error=0.15*r),
        nopt_tot=nopt*tmax
    )

gg <- nopt_new_qmax |> 
    group_by(tmax) |>
    slice_min(qmax, n=1) |>
    ungroup() |>
    pivot_longer(c(nopt_tot_old, nopt_tot)) |> 
    ggplot(aes(qmax, value, col=factor(name, labels=c("full model", "basic model")))) +
    geom_line(linewidth=1.2) +
    scale_color_manual(values=c("darkorange1", "#35A0ABFF")) +
    scale_y_continuous(transform="log", breaks=c(100, 1000, 10000, 50000)) +
    scale_x_continuous(transform="log", breaks=c(0.01, 0.1, 0.3, 0.7), labels=scales::label_percent()) +
    labs(x="prevalence threshold", y="total sequences", col=NULL) +
    theme_half_open() +
    theme(
        legend.position="inside",
        legend.position.inside=c(0.61, 0.84)
    )

ggsave("figs/nopt-full-compare-qmax-1.pdf", width=3.9, height=3.8)

n_max_state <- 214791 / (14*50)

n_max_state <- nopt_new_qmax |> 
    group_by(tmax) |>
    slice_min(qmax, n=1) |>
    ungroup() |>
    mutate(nmax=tmax*n_max_state)

gg +
    geom_line(
        aes(qmax, nmax), data=n_max_state, 
        linewidth=1.2, linetype="dashed", col="gray60",
        inherit.aes=FALSE
    ) +
    theme(legend.position="none")

ggsave("figs/nopt-full-compare-qmax-2.pdf", width=3.9, height=3.8)

###
# hlc <- c("#f5b43d", "#ff856d", "#c63df5", "#f53dbc")

set.seed(1211)

n <- 100
q0 <- 6e-4
r <- 0.12
tmax <- 30

ex_dat <- tibble(
    q=growth_model2(tmax, 0, q0, r),
    y=rbinom(tmax+1, n, q),
    t=seq_along(y)-1
)

r_ml <- map_dbl(c(0, 13), \(t0) {
    optimize(
        \(.r) sum(dbinom(ex_dat$y, n, growth_model2(tmax, t0, q0, .r), log=TRUE)), 
        c(-0.5, 0.5), maximum=TRUE
    )$maximum
})

ml_q <- tibble(
    t=0:tmax,
    q1=growth_model2(tmax, 0, q0, r_ml[1]),
    q2=growth_model2(tmax, 13, q0, r_ml[2]),
) |>
    pivot_longer(-t) |>
    filter(value > 0)

gg <- ggplot(ex_dat, aes(t, y)) +
    geom_point(size=1.6) +
    labs(x="days since sequencing began", y="VOC count") +
    # scale_color_manual(values=hlc[1:3]) +
    # scale_y_continuous(sec.axis=sec_axis(trans=~./100, name="proportion")) +
    # scale_x_continuous(breaks=seq(0, 12, 2)) +
    # coord_cartesian(ylim=c(0, 4.5)) +
    theme_half_open() +
    theme(legend.position="none")

ggsave("figs/unknown-t0-1.pdf", gg, width=4.4, height=2.45)

gg +
    geom_line(aes(y=value*100, col=name), ml_q, linewidth=1.1) +
    geom_text(
        aes(x, y, label=label, col=label), 
        tibble(x=c(10.5, 24), y=c(0.5, 1.5), label=paste0("r = ", round(r_ml, 2))),
        size=5.5
    ) +
    scale_color_manual(values=c("darkorange1", "#35A0ABFF", "darkorange1", "#35A0ABFF"))

ggsave("figs/unknown-t0-2.pdf", width=4.4, height=2.45)

### more realistic example, treating t0 as unknown
n <- 100
q0 <- 1e-4
t0 <- 5.5
r <- 0.12
tmax <- 65

ex_dat <- tibble(
    q=growth_model2(tmax, t0, q0, r),
    y=rbinom(tmax+1, n, q),
    t=seq_along(y)-1
)

ex_dat$q

ll_grid <- expand_grid(t0=seq(-12, 17, length.out=150), r=seq(0.075, 0.175, length.out=100)) |> 
    mutate(loglik=map2_dbl(t0, r, ~gm2_log_lik(c(.x, .y), ex_dat$y, n, tmax, q0)))


(im <- info_mat(tmax, n, t0, q0, r))
cov <- solve(im[1:2, 1:2])
mle_approx <- MASS::mvrnorm(10000, c(r, t0), cov)
colnames(mle_approx) <- c("r", "t0")
mle_approx <- as_tibble(mle_approx)

Z * sqrt(cov[1, 1]) / r * 100
Z * sqrt(var_single_param(n, tmax, t0, q0, r)) / r * 100

# scaling the surface better highlights MORE likely points, so selecting points along
# the band is actually more "fair" in the sense that all of these are relatively very likely
ggplot(ll_grid, aes(r, t0)) +
    geom_raster(aes(fill=1.03^loglik), interpolate=TRUE) +
    stat_ellipse(data=mle_approx, col="darkorange1") +
    stat_ellipse(data=mle_approx, , col="darkorange1", level=0.5, linetype="dotted") +
    annotate("point", x=r, y=5.5, size=1.4, shape=8, col="darkorange1") +
    scale_fill_viridis_c(option="mako") +
    labs(x="growth rate", y="start time", fill="scaled\nlikelihood") +
    scale_x_continuous(expand=expansion()) +
    scale_y_continuous(expand=expansion()) +
    theme(legend.text=element_blank())

ggsave("figs/r-t0-log-lik.pdf", width=3.91, height=3)


###
# r <- 0.07 * 7 # work with a weekly example
# n <- 1000
# 
# set.seed(123)
# sim <- simulate_data(12, n, 3, 1e-4, r, 10000)
# ex_dat <- tibble(prop=sim$y[[1]]/n, t=seq_along(prop)-1)
# 
# res_full <- voc_growth_ci_mc(sim, nsim=NULL, r_int=c(0, 1))
# 
# ml_q <- tibble(
#     t=0:12,
#     q1=growth_model(12, 3, 1e-3, 0.84),
#     q2=growth_model(12, 2, 1e-3, 0.76),
#     q3=growth_model(12, 1, 1e-3, 0.64)
# ) |> 
#     pivot_longer(-t)
#     
# ggplot(ex_dat, aes(t, prop)) +
#     geom_point() +
#     # geom_line(aes(y=value, col=name), ml_q) +
#     labs(x="weeks", y="VOC proportion")
# 
# ggplot(res_full$mle, aes(r, t0)) +
#     geom_hex() +
#     annotate("point", x=r, y=3, col="tomato2", shape=8)
