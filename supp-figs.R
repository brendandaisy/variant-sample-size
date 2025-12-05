library(tidyverse)
library(cowplot)
library(geomtextpath)

source("src/growth-model.R")
source("src/nopt-single-param.R")

pal_discrete <- c("#008080", "#70a494", "#b4c8a8", "#f6edbd", "#ebaa6aff", "#de8a5a", "#ca562cff")

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

ggsave("figs/error-vs-n.pdf", width=2.5, height=2.4)

ggsave("figs/error-asymp-vs-mc-2.pdf", width=5, height=3.5)

# total sequences as function of qmax while varying q0----------------------------

q0 <- 10^(-6:-3)
r <- c(0.05, 0.30)
qmax <- seq(0.01, 0.9, length.out=60)

nopt_qmax_q0 <- expand_grid(qmax, r, q0) |>
    rowwise() |>
    mutate(
        tmax=find_tmax(qmax, q0, r),
        nopt=nopt_single_param(tmax, q0, r, error=r*0.15),
        nopt_tot=nopt*tmax
    ) |> 
    group_by(tmax, r, q0) |>
    slice_min(qmax, n=1) |>
    ungroup()

slice_min(nopt_qmax_q0, nopt_tot)

nopt_qmax_q0 |>
    ggplot(aes(qmax, nopt_tot, group=factor(q0), col=factor(q0))) +
    geom_line(linewidth=1.2) +
    facet_wrap(~r, nrow=1, labeller=label_both) +
    scale_color_manual(values=pal_discrete[c(1, 3, 5, 7)], labels=c("1e-6", "1e-5", "1e-4", "1e-3")) +
    scale_y_continuous(transform="log", breaks=10^(1:4), limits=c(10, NA)) +
    scale_x_continuous(transform="log", breaks=c(0.01, 0.03, 0.10, 0.3, 0.9), labels=scales::label_percent()) +
    # scale_x_continuous(transform="logit", breaks=c(0.01, seq(0.1, 0.9, 0.2), 0.99)) +
    labs(x="prevalence threshold", y="total sequences", col="initial\nprevalence") +
    theme_half_open() +
    theme(
        legend.position="right",
        panel.grid.major=element_line(color="gray70", linetype="dashed"),
        axis.text=element_text(size=rel(0.8))
    )

ggsave("figs/nopt-qmax-varq0.pdf", width=6.5, height=3.6)

# information provided per observation and other analytic derivations-------------
# currently not referenced in the manuscript--------------------------------------
nopt_pct_decrease <- function(t, r, q0) {
    1 - 1 / (1 + info_per_t(t, r, q0) / sum(info_per_t(0:(t-1), r, q0)))
}

q0 <- 1e-3
r <- c(0.07, 0.12)
tmax <- 140

max_info_res <- expand_grid(t=1:tmax, r) |> 
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
    q <- growth_model(tmax, q0, r)
    tibble(r, t=0:tmax, q)
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


#  Prevalence uncertainty and sample size for increasing t------------------------
# r <- 0.12 # daily sequencing for this example (cleaner curve)
# q0 <- 5e-3
# tmax <- find_tmax(0.99, q0, r)
# n <- 6
# 
# var_prev_t <- tibble(t=1:tmax, q=growth_model(tmax, q0, r)[-1]) |>
#     rowwise() |>
#     mutate(
#         vr=var_single_param(t, q0, r)/n,
#         vprev=vr * d_prev_single_param(t, q0, r)^2,
#         ymin=qnorm(0.025, q, sqrt(vprev)), 
#         ymax=qnorm(0.975, q, sqrt(vprev)),
#         dprev=d_prev_single_param(t, q0, r)
#     ) |> 
#     ungroup()
# 
# p1 <- ggplot(var_prev_t) +
#     stat_function(
#         fun=\(x) dnorm(x, r, sqrt(var_single_param(30, q0, r)/n)), 
#         col=pal_discrete[1], linewidth=1.2
#         # geom="area", fill=pal_discrete[1]
#     ) +
#     stat_function(
#         fun=\(x) dnorm(x, r, sqrt(var_single_param(55, q0, r)/n)), 
#         col=pal_discrete[6], linewidth=1.2
#         # geom="area", fill=pal_discrete[1]
#     ) +
#     xlim(r+c(-0.04, 0.04)) +
#     labs(x="growth advantage", y="density") +
#     theme_half_open()
# 
# mdprev <- max(var_prev_t$dprev)
# 
# p2 <- ggplot(var_prev_t, aes(t)) +
#     geom_vline(xintercept=30, linetype="dashed", linewidth=1.2, col=pal_discrete[1]) +
#     geom_vline(xintercept=55, linetype="dashed", linewidth=1.2, col=pal_discrete[6]) +
#     geom_line(aes(y=dprev/mdprev), col="gray50", linewidth=1.1, linetype="dotted") +
#     geom_ribbon(aes(ymin=ymin, ymax=ymax), fill=pal_discrete[3], col=NA, alpha=0.65) +
#     geom_line(aes(y=q), col=pal_discrete[3], linewidth=1.2) +
#     scale_y_continuous(sec.axis=sec_axis(transform=~.*mdprev, name="derivative")) +
#     # scale_x_continuous(breaks=c(0, 5, 10)) +
#     labs(x="days", y="prevalence") +
#     coord_cartesian(ylim=c(0, 1)) +
#     theme_half_open()
# 
# plot_grid(p1, p2, nrow=1, rel_widths=c(0.4, 0.6), labels=c("A", "B"))
# 
# ggsave("figs/delta-method-ex.pdf", width=8, height=3.4)

# Is learning about growth advantage or prevalence more difficult?----------------
# Test the difference for a range of r and q0 settings----------------------------
q0 <- 1e-4
r <- seq(0.05, 0.3, length.out=100)
qmax <- seq(0.01, 0.9999, length.out=100)

nopt_diff_rq <- expand_grid(qmax, r) |> 
    rowwise() |> 
    mutate(
        tmax=find_tmax(qmax, q0, r),
        nopt_prev_per=nopt_g_single_param(tmax, q0, r, d_prev_single_param, error=qmax*0.15),
        nopt_r=nopt_single_param(tmax, q0, r, error=r*0.15),
        diff_tot=tmax*(nopt_prev_per-nopt_r)
    ) |> 
    ungroup() |> 
    mutate(
        diff_tot_scaled=ifelse(diff_tot >= 0, log(diff_tot, 1.5), diff_tot)
    )

ggplot(nopt_diff_rq, aes(qmax, r, fill=diff_tot_scaled)) +
    geom_raster(interpolate=TRUE) +
    scale_fill_gradient2(
        low=pal_discrete[1], high=pal_discrete[6],
        labels=~ifelse(. >= 0, floor(1.5^.), .)
        # transform="log", breaks=scales::breaks_log(3), labels=scales::label_comma()
    ) +
    scale_x_continuous(expand=expansion()) +
    scale_y_continuous(expand=expansion()) +
    labs(
        x="prevalence",
        y="growth advantage",
        fill="total\ndifference"
    )

ggsave("figs/nopt-diff-r-prev.pdf", width=4.1, height=3.2)

nopt_diff_rq |> 
    slice_min(r) |> 
    ggplot(aes(qmax, diff_tot, col=diff_tot_scaled)) +
    geom_line(linewidth=1.15) +
    scale_color_gradient2(
        low=pal_discrete[1], mid="gray90", high=pal_discrete[6],
        labels=~ifelse(. >= 0, floor(1.5^.), .)
        # transform="log", breaks=scales::breaks_log(3), labels=scales::label_comma()
    ) +
    labs(x="prevalence", y="total difference", col="total\ndifference") +
    theme_half_open() +
    theme(panel.grid.major=element_line(color="gray70", linetype="dashed"))

ggsave("figs/nopt-diff-r-prev-2.pdf", width=4.6, height=3.5)

# Sample sizes for difference prevalence and growth combinations------------------
# TODO now likely combined with timing below, and should be cleaned / moved to main
r <- seq(0.05, 0.2, length.out=60)
q0 <- 1e-4
qmax <- c(seq(0.05, 0.95, length.out=60), 0.5)

nopt_prev_grid <- expand_grid(r, qmax) |> 
    rowwise() |> 
    mutate(
        tmax=find_tmax(qmax, q0, r),
        nopt_prev=nopt_g_single_param(tmax, q0, r, d_prev_single_param, error=qmax*0.15),
    ) |> 
    ungroup()

p1 <- nopt_prev_grid |> 
    filter(qmax != 0.5) |> 
    ggplot(aes(qmax, r)) +
    geom_raster(aes(fill=pmax(nopt_prev, 1/7)), interpolate=TRUE) +
    geom_labelvline(
        aes(xintercept=xintercept, label=label, vjust=vjust),
        tibble(xintercept=c(0.05, 0.5, 0.95), label=c("5%", "50%", "95%"), vjust=c(1.3, 0.6, -0.1)), 
        col="gray70", linewidth=0.85
    ) +
    # geom_contour(aes(z=log10(pmax(nopt_prev, 1/7))), col="gray80") +
    scale_fill_viridis_c(
        option="mako", transform="log", 
        breaks=scales::breaks_log(), labels=scales::label_number(1, big.mark=",")
    ) +
    scale_x_continuous(expand=expansion()) +
    scale_y_continuous(expand=expansion()) +
    labs(
        x="target prevalence",
        y="growth advantage",
        fill="sequences\nper day"
    ) +
    theme_half_open() +
    theme(axis.line=element_blank())

np_curves <- nopt_prev_grid |> 
    filter(qmax %in% c(0.05, 0.5, 0.95)) |> 
    mutate(label=str_c(100*qmax, "%"))

pull(np_curves, r)

p2 <- ggplot(np_curves, aes(r, nopt_prev)) +
    # geom_labelpath(aes(label=label), linewidth=1.5, hjust=0.2, straight=TRUE) +
    stat_smooth(
        aes(col=after_stat(y), group=as.factor(qmax), label=label), 
        method="loess", geom="labelpath", linewidth=1.5
    ) +
    scale_y_continuous(
        transform="log", breaks=c(1, 10, 100, 1000), labels=scales::label_comma()
    ) +
    scale_color_viridis_c(option="mako", transform="log") +
    labs(x="growth advantage",y="sequences per day") +
    theme_half_open() +
    theme(legend.position="none")
    
p3 <- ggplot(np_curves, aes(r, nopt_prev*tmax, group=qmax)) +
    geom_labelpath(aes(label=label), linewidth=1.5, straight=TRUE, stat="smooth") +
    scale_y_continuous(transform="log", breaks=scales::breaks_log()) +
    labs(x="growth advantage", y="total sequences") +
    theme_half_open()
    # theme_half_open(10) +
    # theme(plot.background=element_rect(fill="white"))

p2 <- ggdraw(pmain) + draw_plot(pinset, 0.37, 0.45, 0.58, 0.52)

plot_grid(p1, p2, rel_widths=c(1, 0.78), labels="AUTO")

# p2 <- ggplot(nopt_prev_grid, aes(qmax, r)) +
#     geom_raster(aes(fill=tmax*nopt_prev), interpolate=TRUE) +
#     # geom_contour(aes(z=log(nopt_prev)), col="gray80", bins=8) +
#     scale_fill_viridis_c(
#         option="mako", transform="log", 
#         breaks=scales::breaks_log(), labels=scales::label_number(1)
#     ) +
#     scale_x_continuous(expand=expansion()) +
#     scale_y_continuous(expand=expansion()) +
#     labs(
#         x="target prevalence",
#         y="growth advantage",
#         fill="total\nsequences"
#     ) +
#     theme_half_open() +
#     theme(axis.line=element_blank())

plot_grid(p1, p2, labels="AUTO")
ggsave("figs/nopt-prev-grid.pdf", width=7.9, height=3)

# Exploring optimal sample sizes for the time until a prevalence reached----------

## how many sequences to get within 15% for different times and r?
# library(MBA)

q0 <- 1e-4
r <- seq(0.05, 0.20, length.out=200)
qtimes <- c(seq(0.05, 0.95, length.out=100), 0.5)

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

# nt_surf <- mba.surf(select(nopt_timing, time, r, nopt_time), 200, 200, extend=FALSE)
# dimnames(nt_surf$xyz.est$z) <- list(nt_surf$xyz.est$x, nt_surf$xyz.est$y)
# 
# nt_surf <- nt_surf$xyz.est$z |> 
#     as_tibble(rownames="time") |> 
#     pivot_longer(-time) |> 
#     mutate(time=parse_number(time), r=parse_number(name), value=pmax(1/7, value))

nt_path <- nopt_timing |> 
    filter(qtimes %in% c(0.05, 0.5, 0.95)) |> 
    mutate(label=str_c(100*qtimes, "%"))

p4 <- ggplot(nopt_timing, aes(time, r)) +
    geom_line(aes(col=nopt_time, group=r), linewidth=1.2) +
    geom_labelpath(
        aes(label=label), data=nt_path, 
        col="gray70", linewidth=0.85, straight=TRUE
    ) +
    # geom_contour(aes(z=log10(nopt_time)), col="gray80", bins=8) +
    scale_color_viridis_c(
        option="mako", transform="log", 
        breaks=scales::breaks_log(), labels=scales::label_number(1),
        na.value="transparent"
    ) +
    scale_x_continuous(expand=expansion()) +
    scale_y_continuous(expand=expansion()) +
    labs(
        x="target prevalence timing",
        y="growth advantage",
        color="sequences\nper day"
    ) +
    theme_half_open() +
    theme(
        panel.grid.major=element_line(color="gray70", linetype="dashed")
    )

nt_dom <- nopt_timing |> 
    filter(qtimes == 0.5)

p5 <- ggplot(nt_dom, aes(r, nopt_time, col=nopt_time)) +
    geom_line(linewidth=1.5) +
    geom_labelline(data=filter(nt_dom, r > 0.065, r < 0.08), label="Any%", linewidth=1.5, straight=TRUE) +
    # geom_label(x=0.06, y=25, label="Any%", angle=-50) +
    scale_color_viridis_c(option="mako", transform="log") +
    labs(x="growth advantage", y="sequences per day") +
    theme_half_open() +
    theme(
        legend.position="none",
        # panel.grid.major=element_line(color="gray70", linetype="dashed")
    )

p6 <- ggplot(nt_dom, aes(r, nopt_time*time)) +
    geom_labelline(label="Any%", linewidth=1.5, straight=TRUE) +
    scale_y_continuous(transform="log", breaks=scales::breaks_log()) +
    labs(x="growth advantage", y="total sequences") +
    theme_half_open()
    # theme_half_open(10) +
    # theme(plot.background=element_rect(fill="white"))

p4 <- ggdraw(pmain) + draw_plot(pinset, 0.37, 0.45, 0.58, 0.52)

plot_grid(p1, p2, p3, p4, nrow=2, rel_widths=c(1, 0.78), labels="AUTO")

plot_grid(
    p1, p2, p3, p4, p5, p6, nrow=2, ncol=3, 
    rel_widths=c(1, 0.78, 0.78, 1, 0.78, 0.78), labels="AUTO",
    align="v", axis="l"
)

ggsave("figs/nopt-other-var-r-2.pdf", width=12, height=7)

# ggsave("figs/nopt-timing.pdf", width=8.7, height=3.7)
