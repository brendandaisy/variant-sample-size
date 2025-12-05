library(tidyverse)
library(cowplot)

source("src/growth-model.R")
source("src/nopt-single-param.R")

pal_discrete <- c("#008080", "#70a494", "#b4c8a8", "#f6edbd", "#ebaa6aff", "#de8a5a", "#ca562cff")

# Top-line rules of thumb----------------------------------------------------------

# growth advantage
r <- 0.05*7
q0 <- 0.01

nopt_single_param(4, q0, r)
nopt_single_param(4, q0, r, error=0.25*r, level=0.2)
nopt_single_param(8, q0, r)

# prevalence
r <- 0.*7
q0 <- 0.01

(q4wk <- final_prevalence(4, q0, r))
(q8wk <- final_prevalence(8, q0, r))

nopt_g_single_param(4, q0, r, d_prev_single_param, error=0.15*q4wk)
nopt_g_single_param(4, q0, r, d_prev_single_param, error=0.25*q4wk, level=0.2)
nopt_g_single_param(8, q0, r, d_prev_single_param, error=0.15*q8wk)

nopt_single_param(4, q0, r)
nopt_g_single_param(4, q0, r, )
nopt_prev_classic(4, q0, r, error=0.15*final_prevalence(4, q0, r))

nopt_single_param(8, q0, r)
nopt_g_single_param(8, q0, r, d_prev_single_param, error=0.15*final_prevalence(8, q0, r))
nopt_prev_classic(8, q0, r, error=0.15*final_prevalence(8, q0, r))
    
# plots for desired prevalence on log-log scale-----------------------------------
q0 <- 1e-4
# qmax <- c(0.02, 0.1)
# r <- seq(0.05, 0.3, length.out=60)
qmax <- seq(0.01, 0.9, length.out=60)
r <- c(0.01, 0.07, 0.12, 0.18, 0.24, 0.3)

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

p1 <- nopt_grid_qmax |>
    # mutate(label=paste0("q[T]=", qmax)) |>
    filter(nopt >= 1/7) |>
    ggplot(aes(qmax, nopt, group=factor(r), col=factor(r))) +
    geom_line(linewidth=1.2) +
    scale_color_manual(values=pal_discrete[1:6]) +
    scale_y_continuous(transform="log", breaks=c(1, 10, 100, 200)) +
    scale_x_continuous(transform="log", breaks=c(0.01, 0.03, 0.10, 0.3, 0.9), labels=scales::label_percent()) +
    # scale_x_continuous(transform="logit", breaks=c(0.01, seq(0.1, 0.9, 0.2), 0.99)) +
    labs(x="prevalence threshold", y="sequences per day", col="growth\nrate") +
    theme_half_open() +
    theme(
        legend.position="none",
        panel.grid.major=element_line(color="gray70", linetype="dashed"),
        axis.text=element_text(size=rel(0.8))
    )

p2 <- ggplot(nopt_grid_qmax, aes(qmax, nopt_tot, group=factor(r), col=factor(r))) +
    geom_line(linewidth=1.2) +
    # geom_segment(x=log(0.02), xend=log(0.2), y=log(700), yend=log(700) - 1.151*(log(0.2)-log(0.02)), linewidth=1.1, linetype="dotted", col="gray50") +
    # annotate("text", 0.18, 550, label="∇ ≈ -1.15", size=8, col="gray50") +
    scale_color_manual(values=pal_discrete[1:6]) +
    scale_y_continuous(transform="log", breaks=c(1, 100, 500, c(1, 2, 4)*1000)) +
    scale_x_continuous(transform="log", breaks=c(0.01, 0.03, 0.10, 0.3, 0.9), labels=scales::label_percent()) +
    labs(x="prevalence threshold", y="total sequences", col="growth\nadvantage") +
    theme_half_open() +
    theme(
        legend.position="right",
        panel.grid.major=element_line(color="gray70", linetype="dashed"),
        axis.text=element_text(size=rel(0.8))
    )

plot_grid(p1, p2, nrow=1, rel_widths=c(1, 1.33), labels="AUTO")
ggsave("figs/nopt-single-param-qmax.pdf", width=8.1, height=3.6)

## for the reported examples of 3% and 30% bounds in the manuscript:
r_very_small <- 0.0001
q0 <- 1e-4
tm1 <- find_tmax(0.03, q0, r_very_small)
tm2 <- find_tmax(0.3, q0, r_very_small)
tm1*nopt_single_param(tm1, q0, r_very_small, error=r_very_small*0.15)
tm2*nopt_single_param(tm2, q0, r_very_small, error=r_very_small*0.15)

# is adding a given number of sequences better as another day or as spread--------
# over the current interval?------------------------------------------------------

q0 <- 1e-4
qmax_tall <- seq(0.01, 0.99, length.out=100)
r <- c(0.07, 0.12, 0.18, 0.24, 0.3)
# r <- c(0.05, 0.12, 0.3)
# q0 <- c(1e-5, 1e-3)

# pct difference in error when keeping n total constant
err_wait_comp <- expand_grid(n_wide, qmax_tall, r, q0) |> 
    rowwise() |> 
    mutate(
        tmax_tall=find_tmax(qmax_tall, q0, r),
        n_tall=n_wide+1*n_wide/tmax_tall,
        tmax_wide=tmax_tall+1,
        err_tall=qnorm(1 - 0.05/2)*sqrt(var_single_param(tmax_tall, q0, r)/n_tall),
        err_wide=qnorm(1 - 0.05/2)*sqrt(var_single_param(tmax_wide, q0, r)/n_wide),
        pct_diff=(err_wide-err_tall)/err_tall
    ) |> 
    group_by(tmax_tall, r, q0) |>
    slice_min(qmax_tall, n=1) |>
    ungroup()

# rephrase: how many more sequences are necessary with tall strat to achieve the same error?
n_wait_comp <- expand_grid(qmax_tall, r, q0) |> 
    rowwise() |> 
    mutate(
        tmax_tall=find_tmax(qmax_tall, q0, r),
        tmax_wide=tmax_tall+1,
        ntot_tall=nopt_single_param(tmax_tall, q0, r)*tmax_tall,
        ntot_wide=nopt_single_param(tmax_wide, q0, r)*tmax_wide,
        n_diff_tall=ntot_tall-ntot_wide
    ) |> 
    group_by(tmax_tall, r, q0) |>
    slice_min(qmax_tall, n=1) |>
    ungroup()

p1 <- ggplot(err_wait_comp, aes(tmax_tall, pct_diff)) +
    geom_line(aes(col=as.factor(r), group=as.factor(r)), linewidth=1.2) +
    # geom_line(aes(col=r), filter(err_wait_comp, qmax_tall == 0.01), linewidth=1) +
    # geom_vline(
    #     aes(xintercept=max_info_prev), max_info, 
    #     linetype="dashed", col="gray70", linewidth=1.01
    # ) +
    geom_vline(xintercept=0.006, col="white", linewidth=1) +
    geom_vline(xintercept=0, col="gray70", linetype="dashed") +
    # scale_x_continuous(labels=scales::label_percent()) +
    scale_y_continuous(labels=scales::label_percent()) +
    scale_color_manual(values=pal_discrete[2:7]) +
    # scale_color_gradientn(colors=pal_discrete, breaks=c(0.1, 0.2, 0.3)) +
    labs(x="current time", y="error difference", color="growth\nadvantage") +
    theme_half_open() +
    theme(panel.grid.major=element_line(color="gray70", linetype="dashed"), legend.position="none")

p2 <- ggplot(n_wait_comp, aes(tmax_tall, n_diff_tall)) +
    geom_line(aes(col=as.factor(r), group=as.factor(r)), linewidth=1.2) +
    geom_vline(xintercept=0.006, col="white", linewidth=1) +
    geom_vline(xintercept=0, col="gray70", linetype="dashed") +
    # scale_x_continuous(labels=scales::label_percent()) +
    # scale_y_continuous(breaks=seq(0, 8, 2)) +
    scale_color_manual(values=pal_discrete[2:7]) +
    labs(x="current time", y="addtl. total sequences", color="growth\nadvantage") +
    theme_half_open() +
    theme(panel.grid.major=element_line(color="gray70", linetype="dashed"))

plot_grid(p1, p2, nrow=1, rel_widths=c(0.85, 1), labels="AUTO")
ggsave("figs/err-wait-compare-5.pdf", width=8.3, height=3.5)

# retry: only add one day, where the n (spread or added) is increasing
# choice of particular n doesn't matter for % difference
q0 <- 1e-4
qmax_tall <- seq(0.01, 0.99, length.out=60)
r <- seq(0.01, 0.3, length.out=60)
# r <- c(0.05, 0.12, 0.3)
# q0 <- c(1e-5, 1e-3)

# pct difference in error when keeping n total constant
err_wait_comp <- expand_grid(n_wide, qmax_tall, r, q0) |> 
    rowwise() |> 
    mutate(
        tmax_tall=find_tmax(qmax_tall, q0, r),
        n_tall=n_wide+1*n_wide/tmax_tall,
        tmax_wide=tmax_tall+1,
        err_tall=qnorm(1 - 0.05/2)*sqrt(var_single_param(tmax_tall, q0, r)/n_tall),
        err_wide=qnorm(1 - 0.05/2)*sqrt(var_single_param(tmax_wide, q0, r)/n_wide),
        pct_diff=(err_wide-err_tall)/err_tall
    ) |> 
    group_by(tmax_tall, r, q0) |>
    slice_min(qmax_tall, n=1) |>
    ungroup()

# rephrase: how many more sequences are necessary with tall strat to achieve the same error?
n_wait_comp <- expand_grid(qmax_tall, r, q0) |> 
    rowwise() |> 
    mutate(
        tmax_tall=find_tmax(qmax_tall, q0, r),
        tmax_wide=tmax_tall+1,
        ntot_tall=nopt_single_param(tmax_tall, q0, r)*tmax_tall,
        ntot_wide=nopt_single_param(tmax_wide, q0, r)*tmax_wide,
        n_diff_tall=ntot_tall-ntot_wide
        # var_tall=var_single_param(tmax_tall, q0, r),
        # var_wide=var_single_param(tmax_wide, q0, r),
        # n_diff_tall=var_tall/var_wide*tmax_tall - tmax_tall - 1
    ) |> 
    group_by(tmax_tall, r, q0) |>
    slice_min(qmax_tall, n=1) |>
    ungroup()

p1 <- ggplot(err_wait_comp, aes(qmax_tall, pct_diff)) +
    geom_line(aes(col=r, group=r), linewidth=1.6, lineend="round") +
    # geom_line(aes(col=r), filter(err_wait_comp, qmax_tall == 0.01), linewidth=1) +
    # geom_vline(
    #     aes(xintercept=max_info_prev), max_info, 
    #     linetype="dashed", col="gray70", linewidth=1.01
    # ) +
    geom_vline(xintercept=0.006, col="white", linewidth=1) +
    geom_vline(xintercept=0, col="gray70", linetype="dashed") +
    scale_x_continuous(labels=scales::label_percent()) +
    scale_y_continuous(labels=scales::label_percent()) +
    scale_color_gradientn(colors=pal_discrete, breaks=c(0.01, 0.1, 0.2, 0.3)) +
    labs(x="current prevalence", y="error difference", color="growth\nadvantage") +
    theme_half_open() +
    theme(panel.grid.major=element_line(color="gray70", linetype="dashed"), legend.position="none")

p2 <- ggplot(n_wait_comp, aes(qmax_tall, n_diff_tall)) +
    geom_line(aes(col=r, group=r), linewidth=0.85) +
    geom_vline(xintercept=0.006, col="white", linewidth=1) +
    geom_vline(xintercept=0, col="gray70", linetype="dashed") +
    scale_x_continuous(labels=scales::label_percent()) +
    # scale_y_continuous(transform="log", labels=scales::label_log()) +
    scale_color_gradientn(colors=pal_discrete, breaks=c(0.01, 0.1, 0.2, 0.3)) +
    labs(x="current prevalence", y="addtl. total sequences", color="growth\nadvantage") +
    coord_cartesian(ylim=c(0, 250)) +
    theme_half_open() +
    theme(panel.grid.major=element_line(color="gray70", linetype="dashed"))

plot_grid(p1, p2, nrow=1, rel_widths=c(0.85, 1), labels="AUTO")
ggsave("figs/err-wait-compare-4.pdf", width=8.3, height=3.5)

# err_wait_comp |> 
#     ggplot(aes(qmax_tall, pct_diff, col=as.factor(r))) +
#     # geom_vline(aes(xintercept=dom_time), linetype="dashed", col="gray70", linewidth=1.01) +
#     geom_vline(
#         aes(xintercept=max_info_prev), max_info, 
#         linetype="dashed", col="gray70", linewidth=1.01
#     ) +
#     geom_hline(yintercept=0, col="gray70") +
#     geom_line(linewidth=1.1) +
#     facet_grid(r~q0, labeller=label_both) +
#     scale_y_continuous(labels=scales::label_percent()) +
#     scale_color_manual(values=pal_discrete[c(1, 3, 6)]) +
#     labs(x="current days observation", y="error difference", color="growth\nadvantage") +
#     theme_half_open()

# Figure ??: prevalence under periodic sampling------------------------------------
r <- 0.12
q0 <- 5e-3
tmax <- find_tmax(0.99, q0, r)
n <- 6

var_prev_t <- tibble(t=1:tmax, q=growth_model(tmax, q0, r)[-1]) |>
    rowwise() |>
    mutate(
        vr=var_single_param(t, q0, r)/n,
        vprev=vr * d_prev_single_param(t, q0, r)^2,
        ymin=qnorm(0.025, q, sqrt(vprev)), 
        ymax=qnorm(0.975, q, sqrt(vprev)),
        ymin_classic=qnorm(0.025, q, sqrt(q*(1-q)/n)),
        ymax_classic=qnorm(0.975, q, sqrt(q*(1-q)/n)),
        dprev=d_prev_single_param(t, q0, r)
    ) |> 
    ungroup()

pinset <- ggplot(var_prev_t) +
    stat_function(
        fun=\(x) dnorm(x, r, sqrt(var_single_param(30, q0, r)/n)), 
        col=pal_discrete[5], 
        # col="gray60",
        linetype="dotted", linewidth=0.9
        # geom="area", fill=pal_discrete[1]
    ) +
    stat_function(
        fun=\(x) dnorm(x, r, sqrt(var_single_param(55, q0, r)/n)), 
        col=pal_discrete[5], 
        linetype="dashed", linewidth=0.9
        # geom="area", fill=pal_discrete[1]
    ) +
    scale_x_continuous(breaks=c(0.08, 0.12, 0.16), limits=r+c(-0.04, 0.04)) +
    labs(x="est. growth advantage", y="density") +
    theme_half_open(10) +
    theme(
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_rect(fill="white")
    )

mdprev <- max(var_prev_t$dprev)

pmain <- ggplot(var_prev_t, aes(t)) +
    # geom_vline(xintercept=30, linetype="dashed", linewidth=1.1, col=pal_discrete[3]) +
    # geom_vline(xintercept=55, linetype="dashed", linewidth=1.1, col=pal_discrete[5]) +
    geom_segment(x=30, xend=30, y=-Inf, yend=final_prevalence(30, q0, r), linetype="dotted", linewidth=1, col=pal_discrete[5]) +
    geom_segment(x=55, xend=55, y=-Inf, yend=final_prevalence(55, q0, r), linetype="dashed", linewidth=1, col=pal_discrete[5]) +
    geom_ribbon(aes(ymin=ymin_classic, ymax=ymin), fill=pal_discrete[1], col=NA, alpha=0.65) +
    geom_ribbon(aes(ymin=ymax, ymax=ymax_classic), fill=pal_discrete[1], col=NA, alpha=0.65) +
    geom_ribbon(aes(ymin=ymin, ymax=ymax), fill=pal_discrete[7], col=NA, alpha=0.65) +
    geom_line(aes(y=q), col=pal_discrete[7], linewidth=1.2) +
    # scale_y_continuous(sec.axis=sec_axis(transform=~.*mdprev, name="derivative")) +
    scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.0)) +
    labs(x="days", y="estimated prevalence") +
    # coord_cartesian(ylim=c(0, 1)) +
    theme_half_open() +
    theme(plot.margin=unit(c(1, 4, 1, 1), "mm"))

# plot_grid(p1, p2, nrow=1, rel_widths=c(0.4, 0.6), labels=c("A", "B"))
# ggsave("figs/delta-method-ex.pdf", width=4.1, height=3.7)

## Now compare to both cross sectional and the growth advantage
r <- 0.12
q0 <- 1e-4
qmax <- seq(0.01, 0.95, length.out=60)

nopt_prev_qmax <- tibble(qmax) |> 
    rowwise() |> 
    mutate(
        tmax=find_tmax(qmax, q0, r),
        nopt_prev_xs=nopt_prev_classic(tmax, q0, r, error=qmax*0.15),
        nopt_prev_per=nopt_g_single_param(tmax, q0, r, d_prev_single_param, error=qmax*0.15),
        nopt_r=nopt_single_param(tmax, q0, r, error=r*0.15)
    ) |> 
    ungroup()

p2 <- nopt_prev_qmax |> 
    group_by(tmax) |>
    slice_min(qmax, n=1) |>
    ungroup() |>
    pivot_longer(c(nopt_prev_xs, nopt_prev_per, nopt_r)) |> 
    # pivot_longer(contains("nopt")) |> 
    filter(value >= 1/7) |> 
    ggplot(aes(qmax, value, col=name)) +
    geom_line(linewidth=1.12) +
    # geom_vline(xintercept=0.01, linetype="dashed", linewidth=1.04, col="gray50") +
    scale_y_continuous(transform="log", breaks=10^(0:5), labels=scales::label_comma()) +
    scale_x_continuous(
        breaks=c(0.01, 0.25, 0.5, 0.75, 0.99), 
        labels=scales::label_percent()
    ) +
    scale_color_manual(
        values=c(pal_discrete[c(7, 1, 5)]),
        labels=c("prevalence - periodic", "prevalence - cross sect.", "growth advantage")
    ) +
    # xlim(c(0, 180)) +
    labs(x="target prevalence", y="sequences per day", col=NULL) +
    theme_half_open() +
    theme(
        legend.position=c(0.42, 0.84), axis.text=element_text(size=rel(0.8)), 
        legend.box.background=element_rect(fill="white", color="white"),
        panel.grid.major=element_line(color="gray70", linetype="dashed"),
        legend.text=element_text(size=rel(0.85)), plot.margin=unit(c(1, 18, 1, 1), "mm")
    )

p1 <- ggdraw(pmain) + draw_plot(pinset, 0.175, 0.62, 0.37, 0.37)
plot_grid(p1, p2, nrow=1, labels="AUTO", rel_widths=c(0.9888889, 1.1111111))
ggsave("figs/nopt-prev-compare.pdf", width=8.4, height=4)

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
growth_model_bias <- function(tmax, q0, r, b) {
    q <- growth_model(tmax, q0, r)
    q / (q + (1/b)*(1-q))
}

# d_growth_model_bias <- function(tmax, q0, r, b) {
#     a <- 1/q0 - 1
#     ts <- 0:tmax
#     (a*b*ts*exp(-r*ts))/(1 + a*b*exp(-r*ts))^2
# }
# 
# nopt_single_param_bias <- function(tmax, q0, r, b, error=0.05, level=0.95) {
#     q_bias <- growth_model_bias(tmax, 0, q0, r, b)
#     dq_bias <- d_growth_model_bias(tmax, q0, r, b)
#     It <- dq_bias^2 * (1/(1-q_bias) + 1/q_bias)
#     Z <- -qnorm((1 - level)/2)
#     return((Z/error)^2 / sum(It))
# }
# 
# find_tmax_bias <- function(qmax, t0, q0, r, b) {
#     q_bias <- growth_model_bias(10000, t0, q0, r, b)
#     min(which(q_bias >= qmax)) - 1
# }

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
