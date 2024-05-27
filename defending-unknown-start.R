# analysis to justify our decision to consider the start time unknown
# a good supp fig 1
# 1. a figure illustrating the data generation process, motivating: "when did we start?"
# 2. we could analyse proportional data with 3 models: linear transform, simple logistic, logistic w/ unknown start

# supp fig 2: illustration of patching the logistic with unknown start, with q as func of t0

library(tidyverse)
library(cowplot)

source("src/growth-model.R")

hlc <- c("#f5b43d", "#ff856d", "#c63df5", "#f53dbc")

n <- 100
tmax <- 12

ex_dat <- tibble(
    y=c(0, 0, 0, 0, 0, 0, 2, 0, 1, 2, 2, 4, 1),
    t=seq_along(y)-1
)

ll_grid <- expand_grid(t0=seq(0, 6, 0.1), r=7*seq(-0.1, 0.2, 0.01)) |> 
    mutate(loglik=map2_dbl(t0, r, ~gm2_log_lik(c(.x, .y), ex_dat$y, 100, tmax, 1e-4)))

ml_q <- tibble(
    t=0:tmax,
    q1=growth_model2(tmax, 0.1, 1e-4, 0.52),
    # q2=growth_model2(tmax, 2, 1e-4, 0.61),
    q2=growth_model2(tmax, 3.5, 1e-4, 0.77),
    q3=growth_model2(tmax, 5.7, 1e-4, 1),
) |> 
    pivot_longer(-t) |> 
    filter(value > 0)

points <- tribble(
    ~t0, ~r,
    0.1, 0.52,
    # 2, 0.61,
    3.5, 0.77,
    5.7, 1
)

ggplot(ex_dat, aes(t, y)) +
    geom_line(aes(y=value*100, col=name), ml_q, linewidth=1.1) +
    geom_point(size=1.6) +
    labs(x="week", y="number positive") +
    scale_color_manual(values=hlc[1:3]) +
    scale_y_continuous(sec.axis=sec_axis(trans=~./100, name="proportion")) +
    scale_x_continuous(breaks=seq(0, 12, 2)) +
    coord_cartesian(ylim=c(0, 4.5)) +
    theme_half_open() +
    theme(legend.position="none")

ggsave("figs/example-fits.pdf", width=4.4, height=2.45)

# scaling the surface better highlights MORE likely points, so selecting points along
# the band is actually more "fair" in the sense that all of these are relatively very likely
ggplot(ll_grid, aes(r, t0)) +
    geom_raster(aes(fill=1.02^loglik), interpolate=TRUE) +
    geom_point(data=points, col=hlc[1:3], size=1.8, shape=8) +
    scale_fill_viridis_c(option="mako") +
    labs(x="selection coefficient", y="start time") +
    scale_x_continuous(expand=expansion()) +
    scale_y_continuous(expand=expansion()) +
    theme(legend.position="none")

ggsave("figs/example-log-lik.pdf", width=3, height=3)

r <- 0.07 * 7 # work with a weekly example
n <- 1000

set.seed(123)
sim <- simulate_data(12, n, 3, 1e-4, r, 10000)
ex_dat <- tibble(prop=sim$y[[1]]/n, t=seq_along(prop)-1)

res_full <- voc_growth_ci_mc(sim, nsim=NULL, r_int=c(0, 1))

ml_q <- tibble(
    t=0:12,
    q1=growth_model(12, 3, 1e-3, 0.84),
    q2=growth_model(12, 2, 1e-3, 0.76),
    q3=growth_model(12, 1, 1e-3, 0.64)
) |> 
    pivot_longer(-t)
    
ggplot(ex_dat, aes(t, prop)) +
    geom_point() +
    # geom_line(aes(y=value, col=name), ml_q) +
    labs(x="weeks", y="VOC proportion")

ggplot(res_full$mle, aes(r, t0)) +
    geom_hex() +
    annotate("point", x=r, y=3, col="tomato2", shape=8)
