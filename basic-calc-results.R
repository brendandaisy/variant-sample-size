library(tidyverse)
library(cowplot)

source("src/growth-model.R")
source("src/nopt-single-param.R")

pal_discrete <- c("#008080", "#70a494", "#b4c8a8", "#f6edbd", "#ebaa6aff", "#de8a5a", "#ca562cff")

# Panel A: example logistic curves with different selection coefficients----------
q0 <- 1e-4
# qmax <- c(0.01, 0.1)
r <- c(0.07, 0.12)
# params <- expand_grid(qmax, r)

curves <- map_dfr(r, \(r) {
    q <- growth_model(120, q0, r)
    obs10 <- rbinom(length(q), 10, q) / 10
    obs100 <- rbinom(length(q), 100, q) / 100
    tibble(r, t=0:120, q, `10`=obs10, `100`=obs100)
})

## version for Inkscape
# curves |> 
#     filter(r == 0.12) |> 
#     ggplot(aes(t, q)) +
#     geom_line(col="gray60", linewidth=1.6) +
#     theme_half_open() +
#     theme(
#         axis.title=element_blank(),
#         axis.text=element_blank(),
#         axis.ticks=element_blank(),
#         axis.line=element_line(color="gray70", linewidth=1)
#     )
# 
# ggsave("figs/logistic-growth-blank.pdf", width=4.5, height=2.9)

p_text <- function() theme(
    text=element_text(size=12),
    axis.text=element_text(size=10),
    legend.text=element_text(size=10),
)

p1 <- curves |> 
    pivot_longer(c(`10`, `100`)) |> 
    ggplot(aes(t, q, col=factor(r), group=factor(r))) +
    geom_line() +
    geom_point(aes(y=value, shape=name), alpha=0.6) +
    scale_shape_manual(values=c(1, 16)) +
    scale_color_manual(values=c(pal_discrete[1], pal_discrete[6])) +
    labs(
        x="days since variant introduction", y="variant proportion",
        col="growth\nadvantage", shape="sample\nsize"
    ) +
    guides(color=guide_legend(order=1)) +
    theme_half_open() +
    p_text() +
    theme(,
        legend.key.height=unit(0.8, "cm"), 
        plot.margin=unit(c(10, 0, 5, 0), "mm"),
        legend.key.spacing.y=unit(-2.5, "mm"),
        legend.spacing=unit(1, "mm")
    )

## nopt for different settings
nopt_grid1 <- expand_grid(tmax=10:70, r=seq(0.05, 0.3, length.out=60)) |>
    rowwise() |>
    mutate(
        nopt=max(nopt_single_param(tmax, q0=1e-4, r, error=r*0.15), 1),
        nopt_total=nopt*tmax,
    )

nopt_grid2 <- expand_grid(r=seq(0.05, 0.3, length.out=60), qm=seq(-5, -2, length.out=60)) |>
    rowwise() |>
    mutate(
        nopt=max(nopt_single_param(tmax=30, q0=10^qm, r, error=r*0.15), 1),
        nopt_total=nopt*tmax,
    )

nopt_grid3 <- expand_grid(tmax=10:70, qm=seq(-5, -2, length.out=60)) |>
    rowwise() |>
    mutate(
        nopt=max(nopt_single_param(tmax, q0=10^qm, r=0.12, error=r*0.15), 1),
        nopt_total=nopt*tmax,
    )

pal_lim <- max(max(nopt_grid1$nopt), max(nopt_grid2$nopt), max(nopt_grid3$nopt))

p2 <- ggplot(nopt_grid1, aes(tmax, r)) +
    geom_raster(aes(fill=nopt), interpolate=TRUE) +
    geom_contour(aes(z=log10(nopt)), col="gray80", bins=8) +
    # geom_hline(yintercept=r, col="gray70", linetype="31", linewidth=1.05) +
    scale_fill_viridis_c(option="mako", transform="log", breaks=scales::breaks_log(7), labels=scales::label_comma(), limits=c(NA, pal_lim)) +
    scale_x_continuous(expand=expansion()) +
    scale_y_continuous(expand=expansion()) +
    labs(
        x="days of observation", 
        y="growth advantage", 
        fill="sequences\nper day"
    ) +
    p_text() +
    theme(legend.position="none")

p3 <- ggplot(nopt_grid2, aes(qm, r)) +
    geom_raster(aes(fill=nopt), interpolate=TRUE) +
    geom_contour(aes(z=log10(nopt)), col="gray80", bins=8) +
    # geom_hline(yintercept=r, col="gray70", linetype="31", linewidth=1.05) +
    scale_fill_viridis_c(option="mako", transform="log", breaks=scales::breaks_log(7), labels=scales::label_comma(), limits=c(NA, pal_lim)) +
    scale_x_continuous(expand=expansion(), labels=c("1e-5", "1e-4", "1e-3", "1e-2")) +
    scale_y_continuous(expand=expansion()) +
    labs(
        x="initial proportion", 
        y="growth advantage", 
        fill="sequences per day"
    ) +
    p_text() +
    theme(legend.position="none")

p4 <- ggplot(nopt_grid3, aes(tmax, qm)) +
    geom_raster(aes(fill=nopt), interpolate=TRUE) +
    geom_contour(aes(z=log10(nopt)), col="gray80", bins=8) +
    # geom_hline(yintercept=r, col="gray70", linetype="31", linewidth=1.05) +
    scale_fill_viridis_c(option="mako", transform="log", breaks=scales::breaks_log(7), labels=scales::label_comma()) +
    scale_x_continuous(expand=expansion()) +
    scale_y_continuous(expand=expansion(), labels=c("1e-5", "1e-4", "1e-3", "1e-2")) +
    labs(
        x="days of observation", 
        y="initial proportion", 
        fill="sequences\nper day"
    ) +
    p_text() +
    theme(legend.position="none", plot.background=element_rect(fill="transparent", color=NA))

legend <- get_legend(
    p2 + theme(
        legend.position="right", 
        legend.key.height=unit(0.7, "cm"),
        legend.margin=margin(0, 5, 15, 0, "mm")
    )
)

plot_grid(plot_grid(
    plot_grid(NULL, p1, NULL, nrow=1, rel_widths=c(0.1, 0.8, 0.1), labels=c("", "A", "")),
    plot_grid(
        p2, p3, p4, NULL, legend,
        nrow=1,
        rel_widths=c(1, 1, 1, 0.25, 0.09),
        labels=c("B", "C", "D", "", "")
    ),
    nrow=2, rel_heights=c(1, 0.8)
), NULL, nrow=1, rel_widths=c(1, 0.05))

ggsave("figs/nopt-r-first-plot.pdf", width=7.2, height=5)

nopt_grid2 |>
    filter(isTRUE(all.equal(r, 0.10084746, tol=0.01)) | r == 0.3) |> 
    ggplot(aes(qm, nopt, col=r, group=r)) +
    geom_line()

nopt_grid4 |>
    ggplot(aes(rel_err, nopt, col=r, group=r)) +
    geom_line()

plot_grid(p1, p2, nrow=1, rel_widths=c(0.55, 0.45), labels=c("A", "B"))

ggsave("figs/nopt-first-fig.pdf", width=4.7, height=3.2)

# information provided per observation and other analytic derivations-------------
info_per_t <- function(t, r, q0) {
    a <- 1/q0 - 1
    eterm <- exp(-r*t)
    (t^2*a*eterm) / (1 + a*eterm)^2
}

nopt_pct_decrease <- function(t, r, q0) {
    1 - 1 / (1 + info_per_t(t, r, q0) / sum(info_per_t(0:(t-1), r, q0)))
}

# solve for where the derivative of I_t(r) == 0 in t
# f is said derivative
info_root_t <- function(r, q0) {
    a <- 1/q0 - 1
    tmax <- find_tmax(0.99, q0, r)
    f <- function(t) (a*t*exp(r*t) * (a*(r*t + 2) + exp(r*t)*(2 - r*t)))/(a + exp(r*t))^3
    uniroot(f, c(1, tmax), tol=0.00001)
}

q0 <- 1e-4
r <- c(0.07, 0.12)

max_info_res <- expand_grid(t=1:170, r) |> 
    rowwise() |> 
    mutate(
        dom=dom_time(q0, r),
        info=info_per_t(t, r, q0),
        nopt=nopt_single_param(t, q0, r, error=0.15*r),
        pct_decrease=map_dbl(t, nopt_pct_decrease, r, q0)
    ) |> 
    mutate(r=as.factor(r))

tail(max_info_res)

ymax <- max(max_info_res$info)

p1 <- map_dfr(r, \(r) {
    q <- growth_model(170, q0, r)
    tibble(r, t=0:170, q)
}) |> 
    ggplot(aes(t, q, col=factor(r), group=factor(r))) +
    geom_vline(aes(xintercept=dom, col=r), max_info_res, linetype="dashed") +
    geom_line(linewidth=1.05) +
    scale_color_manual(values=c(pal_discrete[1], pal_discrete[6])) +
    labs(x=NULL, y="variant\nproportion") +
    theme_half_open() +
    theme(legend.position="none")

p2 <- ggplot(max_info_res, aes(t, col=r, group=r)) +
    geom_vline(aes(xintercept=dom, col=r), linetype="dashed", show.legend=FALSE) +
    geom_line(aes(y=info), linewidth=1.05) +
    scale_color_manual(values=c(pal_discrete[1], pal_discrete[6])) +
    labs(
        x=NULL, 
        y="information\nper sequence",
        col="growth\nadvantage"
    ) +
    theme_half_open() +
    theme(legend.position="none")

p3 <- ggplot(max_info_res, aes(t, col=r, group=r)) +
    geom_vline(aes(xintercept=dom, col=r), linetype="dashed", show.legend=FALSE) +
    geom_line(aes(y=nopt), linewidth=1.05) +
    scale_color_manual(values=c(pal_discrete[1], pal_discrete[6])) +
    scale_y_log10(n.breaks=6) +
    # coord_cartesian(ylim=c(0, 100)) +
    labs(
        x=NULL, 
        y="sample\nsize"
    ) +
    theme_half_open() +
    theme(legend.position="none")
    
p4 <- ggplot(max_info_res, aes(t, col=r, group=r)) +
    geom_vline(aes(xintercept=dom, col=r), linetype="dashed", show.legend=FALSE) +
    geom_line(aes(y=pct_decrease), linewidth=1.05) +
    scale_color_manual(values=c(pal_discrete[1], pal_discrete[6])) +
    scale_y_continuous(labels=scales::label_percent()) +
    labs(
        x="days of observation", 
        y="% decrease\nsample size",
        col="growth\nadvantage"
    ) +
    theme_half_open() +
    theme(legend.position="none")

legend <- get_legend(p2 + theme(legend.position="right"))

plot_grid(
    plot_grid(p1, p2, p3, p4, nrow=4, align="v", axis="l", labels="AUTO"),
    legend,
    nrow=1, rel_widths=c(1, 0.2)
)

ggsave("figs/info-per-t.pdf", width=5.7, height=8)
    
# plots for desired prevalence on log-log scale-----------------------------------
q0 <- 1e-4
# qmax <- c(0.02, 0.1)
# r <- seq(0.05, 0.3, length.out=60)
qmax <- seq(0.01, 0.9, length.out=60)
r <- c(0.0001, 0.07, 0.12, 0.18, 0.24, 0.3)

nopt_grid_qmax <- expand_grid(qmax, r) |>
    rowwise() |>
    mutate(
        tmax=find_tmax(qmax, q0, r),
        nopt=nopt_single_param(tmax, q0, r, error=r*0.15),
        nopt_tot=nopt*tmax
    ) |> 
    group_by(tmax, r) |>
    slice_min(qmax, n=1) |>
    ungroup()


lm(log(nopt) ~ log(qmax), filter(nopt_grid_qmax, qmax <= 0.5))

# percentage decrease in daily samples needed for a 1% increase in prevalence
(1.01^(-1.522) - 1) * 100
# percentage decrease in total samples needed for a 1% increase in prevalence
(1.01^(-1.151) - 1) * 100

lim <- expand_grid(qmax, r=1e-6) |>
    rowwise() |>
    mutate(
        tmax=find_tmax(qmax, q0, r),
        nopt=nopt_single_param(tmax, q0, r, error=r*0.15),
        nopt_tot=nopt*tmax
    ) |> 
    group_by(tmax, r) |>
    slice_min(qmax, n=1) |>
    ungroup()

lim

find_tmax(0.01, q0, 1e-6)*nopt_single_param(find_tmax(0.01, q0, 1e-6), q0, 1e-6, error=1e-6*0.15)

p1 <- nopt_grid_qmax |>
    # mutate(label=paste0("q[T]=", qmax)) |>
    filter(nopt >= 1/7) |>
    ggplot(aes(qmax, nopt, group=factor(r), col=factor(r))) +
    geom_line(linewidth=1.2) +
    scale_color_manual(values=pal_discrete[1:6]) +
    scale_y_continuous(transform="log", breaks=c(1, 10, 100, 200)) +
    scale_x_continuous(transform="log", breaks=c(0.01, 0.03, 0.10, 0.5, 0.9), labels=scales::label_percent()) +
    # scale_x_continuous(transform="logit", breaks=c(0.01, seq(0.1, 0.9, 0.2), 0.99)) +
    labs(x="prevalence threshold", y="sequences per day", col="growth\nrate") +
    theme_half_open() +
    theme(
        legend.position="none",
        axis.text=element_text(size=rel(0.8))
    )

p2 <- ggplot(nopt_grid_qmax, aes(qmax, nopt_tot, group=factor(r), col=factor(r))) +
    geom_line(linewidth=1.2) +
    # geom_segment(x=log(0.02), xend=log(0.2), y=log(700), yend=log(700) - 1.151*(log(0.2)-log(0.02)), linewidth=1.1, linetype="dotted", col="gray50") +
    # annotate("text", 0.18, 550, label="∇ ≈ -1.15", size=8, col="gray50") +
    scale_color_manual(values=pal_discrete[1:6]) +
    scale_y_continuous(transform="log", breaks=c(1, 100, 500, c(1, 2, 4)*1000)) +
    scale_x_continuous(transform="log", breaks=c(0.01, 0.03, 0.10, 0.5, 0.9), labels=scales::label_percent()) +
    labs(x="prevalence threshold", y="total sequences", col="growth\nadvantage") +
    theme_half_open() +
    theme(
        legend.position="right",
        axis.text=element_text(size=rel(0.8))
    )

plot_grid(p1, p2, nrow=1, rel_widths=c(1, 1.33), labels="AUTO")
ggsave("figs/nopt-single-param-qmax.pdf", width=8.1, height=3.6)

# key points:
# waiting until VOC has grown will significantly reduce the number of samples needed (per day)
# however, the total number of samples you will need for a particular prevalence threshold is pretty 
# much the same regardless of how fast the variant grows

#  variant proportion over time for single param model----------------------------
r <- 0.12
q0 <- 1e-4
qmax <- seq(0.01, 0.95, length.out=100)

nopt_single_param(60, q0, r, error=r*0.15)
nopt_g_single_param(60, q0, r, d_dom_single_param, error=dom_time(q0, r)*0.15)

nopt_oth_qmax <- tibble(qmax) |> 
    rowwise() |> 
    mutate(
        tmax=find_tmax(qmax, q0, r),
        nopt_prev_xs=nopt_prev_classic(tmax, q0, r, error=qmax*0.15),
        nopt_prev_per=nopt_g_single_param(tmax, q0, r, d_prev_single_param, error=qmax*0.15),
        nopt_r=nopt_single_param(tmax, q0, r, error=r*0.15),
        nopt_dom=nopt_g_single_param(tmax, q0, r, d_dom_single_param, error=dom_time(q0, r)*0.15),
        nopt_concern=nopt_g_single_param(
            tmax, q0, r,
            \(t, q0, r) d_tmax_single_param(0.05, t, q0, r),
            error=find_tmax(0.05, q0, r)*0.15
        )
    ) |> 
    ungroup()

p1 <- nopt_oth_qmax |> 
    group_by(tmax) |>
    slice_min(qmax, n=1) |>
    ungroup() |>
    pivot_longer(c(nopt_prev_xs, nopt_prev_per, nopt_r)) |> 
    # pivot_longer(contains("nopt")) |> 
    filter(value >= 1/7) |> 
    ggplot(aes(qmax, value*tmax, col=name)) +
    geom_line(linewidth=1.06) +
    # geom_vline(xintercept=0.01, linetype="dashed", linewidth=1.04, col="gray50") +
    scale_y_continuous(transform="log", breaks=10^(0:5), labels=scales::label_comma()) +
    scale_x_continuous(
        breaks=c(0.01, 0.25, 0.5, 0.75, 0.99), 
        labels=scales::label_percent()
    ) +
    scale_color_manual(
        values=c(pal_discrete[c(7, 1, 5)]),
        labels=c("prevalence - periodic", "prevalence - cross sect.", "growth advantage/\nwaiting time")
    ) +
    # xlim(c(0, 180)) +
    labs(x="prevalence", y="total sequences needed", col=NULL) +
    theme_half_open() +
    theme(
        legend.position=c(0.42, 0.84), axis.text=element_text(size=rel(0.8)), 
        legend.text=element_text(size=rel(0.85)), plot.margin=unit(c(1, 18, 1, 1), "mm")
    )

### how many sequences to get within 1 day for different times and r?
q0 <- 1e-4
qmax <- 0.05
r <- seq(0.07, 0.15, length.out=60)
qtimes <- seq(0.05, 0.95, length.out=60)

nopt_timing <- expand_grid(r, qtimes) |> 
    rowwise() |> 
    mutate(
        time=find_tmax(qtimes, q0, r),
        nopt_time=nopt_g_single_param(
            tmax, q0, r,
            \(t, q0, r) d_tmax_single_param(qtimes, t, q0, r),
            # error=7
            error=0.15*time
        )
    ) |> 
    ungroup()

library(MBA)

nt_surf <- mba.surf(select(nopt_timing, time, r, nopt_time), 100, 100)
dimnames(nt_surf$xyz.est$z) <- list(nt_surf$xyz.est$x, nt_surf$xyz.est$y)

nt_surf <- nt_surf$xyz.est$z |> 
    as_tibble(rownames="time") |> 
    pivot_longer(-time) |> 
    mutate(time=parse_number(time), r=parse_number(name))

p2 <- ggplot(nt_surf, aes(time, r)) +
    geom_tile(aes(fill=value, col=value)) +
    # geom_contour(aes(z=log10(nopt_time)), col="gray80", bins=8) +
    scale_fill_viridis_c(
        option="mako", transform="log", 
        breaks=scales::breaks_log(), labels=scales::label_comma(),
        na.value="transparent"
    ) +
    scale_color_viridis_c(
        option="mako", transform="log", 
        guide="none",
        na.value="transparent"
    ) +
    scale_y_continuous(expand=expansion()) +
    labs(
        x="time until 5%-95% prevalence reached",
        y="growth advantage",
        fill="sequences\nper day"
    ) +
    theme_half_open() +
    background_grid(major="xy")

plot_grid(p1, p2, nrow=1, align="h", labels="AUTO")
ggsave("figs/nopt-other-metrics-res-v2.pdf", width=8.5, height=3.7)

### how many sequences to recognize the time of X%, by the time the variant is at X%?
# might be something interesting here??

# increases with r because, to reach a relative accuracy of 0.15*time, is same n* as
# to reach relative accuracy of 0.15r. Thus, when r decreases, more samples are required to 
# reach same accuracy

# r <- 0.07
# q0 <- 1e-4
# qmax <- seq(0.01, 0.95, length.out=100)
# 
# nopt_roll_time <- tibble(qmax) |> 
#     rowwise() |> 
#     mutate(
#         time=find_tmax(qmax, q0, r),
#         nopt_timing=nopt_g_single_param(
#             tmax, q0, r,
#             \(t, q0, r) d_tmax_single_param(qmax, t, q0, r),
#             error=time*0.15
#         )
#     ) |> 
#     ungroup()
# 
# nopt_roll_time |> 
#     group_by(time) |>
#     slice_min(qmax, n=1) |>
#     ungroup() |>
#     mutate(nopt_timing_tot=time*nopt_timing) |> 
#     pivot_longer(contains("nopt")) |> 
#     ggplot(aes(qmax)) +
#     geom_line(aes(y=value)) +
#     facet_wrap(~name, nrow=1, scales="free")

#  nopt by prevalance threshold, including sampling bias--------------------------
growth_model_bias <- function(tmax, t0, q0, r, b) {
    q <- growth_model2(tmax, t0, q0, r)
    q / (q + b*(1-q))
}

growth_model_bias2 <- function(tmax, q0, r, b) {
    a <- 1/q0 - 1
    ts <- 0:tmax # this one only works when t0=0 I think
    1 / (1 + a*b*exp(-r * ts))
}

d_growth_model_bias <- function(tmax, q0, r, b) {
    a <- 1/q0 - 1
    ts <- 0:tmax
    (a*b*ts*exp(-r*ts))/(1 + a*b*exp(-r*ts))^2
}

nopt_single_param_bias <- function(tmax, q0, r, b, error=0.05, level=0.95) {
    q_bias <- growth_model_bias(tmax, 0, q0, r, b)
    dq_bias <- d_growth_model_bias(tmax, q0, r, b)
    It <- dq_bias^2 * (1/(1-q_bias) + 1/q_bias)
    Z <- -qnorm((1 - level)/2)
    return((Z/error)^2 / sum(It))
}

find_tmax_bias <- function(qmax, t0, q0, r, b) {
    q_bias <- growth_model_bias(10000, t0, q0, r, b)
    min(which(q_bias >= qmax)) - 1
}

# g <- growth_model(120, 0, q0, 0.12)
# q <- growth_model_bias(120, 0, q0, 0.12, 1/3)
# 
# plot(log(g / (1 - g)))
# points(log(q / (1 - q)), col="blue")

q0 <- 1e-4
qmax <- c(0.01, 0.05, 0.1)
tmax <- c(30, 60)
r <- c(0.07, 0.12)
b <- c(1/3, 1, 3)
# b <- seq(1/3, 3, length.out=100) # over sampled/enrichment, neutral, over sampled

nopt_res_bias <- expand_grid(r, tmax, b) |> 
    rowwise() |> 
    mutate(
        # tmax=find_tmax_bias(qmax, 0, q0, r, b),
        nopt=ceiling(nopt_single_param_bias(tmax, q0, r, b, error=0.25*r)),
        nopt_total=tmax * nopt
    ) |> 
    ungroup()

# for now just inspect a few examples:
nopt_res_bias

nopt_res_bias |> 
    ggplot(aes(1/b, nopt, col=factor(r), group=factor(r))) +
    geom_line() +
    scale_y_sqrt()

nopt_res_bias |> 
    filter(r == 0.07, qmax == 0.05) |> 
    pivot_longer(c(nopt, nopt_total)) |> 
    ggplot(aes(1/b, value)) +
    # geom_tile(aes(fill=nopt_total)) +
    # geom_text(aes(label=nopt_total))
    geom_line() +
    facet_wrap(~name, scales="free_y")
