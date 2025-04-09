library(tidyverse)
library(cowplot)

source("src/growth-model.R")
source("src/nopt-single-param.R")

pal_discrete <- c("#008080", "#70a494", "#b4c8a8", "#f6edbd", "#ebaa6aff", "#de8a5a", "#ca562cff")

# TODO: all of these needs to be edited for new nopt functions (removing t0 and whatnot)

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

gg <- curves |> 
    pivot_longer(c(`10`, `100`)) |> 
    ggplot(aes(t, q, col=factor(r), group=factor(r))) +
    geom_line() +
    geom_point(aes(y=value, shape=name), alpha=0.6) +
    scale_shape_manual(values=c(1, 16)) +
    scale_color_manual(values=c("#ca562c", "#008080")) +
    labs(
        x="days since variant introduction", y="variant proportion",
        col="growth\nrate", shape="sample\nsize"
    ) +
    guides(color=guide_legend(order=1)) +
    theme_half_open()

ggsave("figs/logistic-growth.pdf", gg, width=5.5, height=4.4)

gg +
    theme(
        legend.position="none",
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_rect(fill = "transparent",
                                      colour = NA_character_),
        plot.background=element_rect(fill = "transparent",
                                      colour = NA_character_)
    )

ggsave("figs/logistic-growth-title-page.pdf", width=6, height=4.1, bg="transparent")

## version for Inkscape
curves |> 
    filter(r == 0.12) |> 
    ggplot(aes(t, q)) +
    geom_line(col="gray60", linewidth=1.6) +
    theme_half_open() +
    theme(
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_line(color="gray70", linewidth=1)
    )

ggsave("figs/logistic-growth-blank.pdf", width=4.5, height=2.9)

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

# Panel B: comparing asymptotic and Monte Carlo error as n increases--------------
q0 <- 1e-4
r <- 0.07
t0 <- 0 # just in case..
tmax <- find_tmax(c(0.02, 0.1), 0, q0, r)

Z <- -qnorm((1 - 0.95)/2)
err_form <- expand_grid(tmax, ns=seq(1, 90, 1)) |> 
    rowwise() |> 
    mutate(
        n=ns,
        err_raw={
            vars <- map_dbl(r, ~var_single_param(n, tmax, 0, q0, .))
            Z * sqrt(vars)
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
ggplot(err_form, aes(n, err, col=factor(tmax, labels=c(0.02, 0.1)))) +
    geom_hline(yintercept=25, col="gray70", linewidth=1.2, linetype="dashed") +
    geom_hline(yintercept=15, col="gray70", linewidth=1.2, linetype="dotted") +
    geom_point(aes(y=err_mc), err_mc) +
    geom_line(aes(y=err_mc), err_mc, linetype="dashed") +
    geom_line(linewidth=1.2) +
    scale_y_continuous(labels=scales::label_percent(scale=1)) +
    scale_color_manual(values=c("darkorange1", "#35A0ABFF")) +
    labs(x="sequences per day", y="error", col="final\nprevalence") +
    coord_cartesian(ylim=c(0, 100)) +
    theme_half_open()

ggsave("figs/error-asymp-vs-mc-2.pdf", width=5, height=3.5)

### plot 
q0 <- 1e-4
# n <- 100

vgrid <- expand_grid(tmax=5:70, r=seq(0.05, 0.3, length.out=60)) |>
    rowwise() |>
    mutate(
        nopt=max(nopt_single_param(tmax, q0, r, error=r*0.15), 1),
        nopt_total=nopt*tmax,
    )

ggplot(vgrid, aes(tmax, r)) +
    geom_raster(aes(fill=nopt), interpolate=TRUE) +
    geom_contour(aes(z=log10(nopt)), col="gray80", bins=8) +
    # geom_hline(yintercept=r, col="gray70", linetype="31", linewidth=1.05) +
    scale_fill_viridis_c(option="mako", transform="log10", breaks=scales::breaks_log(7)) +
    scale_x_continuous(expand=expansion()) +
    scale_y_continuous(expand=expansion()) +
    labs(
        x="days since variant introduction", 
        y="growth advantage", 
        fill="sequences per day"
    )

ggsave("figs/nopt-single-param-tmax.pdf", width=5.2, height=3.7)

# information provided per observation and other analytic derivations-------------
info_per_t <- function(t, r, q0) {
    a <- 1/q0 - 1
    eterm <- exp(-r*t)
    (t^2*a*eterm) / (1 + a*eterm)^2
}

nopt_pct_decrease <- function(t, r, q0) {
    1 - 1 / (1 + info_per_t(t, r, q0) / sum(info_per_t(0:(t-1), r, q0)))
}

info_root_t <- function(r, q0) {
    a <- 1/q0 - 1
    tmax <- find_tmax(0.99, 0, q0, r)
    f <- function(t) (a*t*exp(r*t) * (a*(r*t + 2) + exp(r*t)*(2 - r*t)))/(a + exp(r*t))^3
    uniroot(f, c(1, tmax), tol=0.00001)
}

max_decel_time <- function(r, q0) {

}

## notes from 2/14:
## found a very interesting theoretical thing that the most informative (information maxing)
## time point is actually just after the domination time!!
## struggled with deriving this critical point however. 
## overall, come back to this plot as an interesting way to show how much new information an
## extra day of observation adds, and why

q0 <- 1e-4
r <- 0.07

tibble(
    t=0:170, info=info_per_t(t, r, q0), 
    pct_decrease=map_dbl(t, nopt_pct_decrease, r, q0)
) |> 
    ggplot(aes(t, info)) +
    geom_vline(xintercept=dom_time(0, q0, r)) +
    geom_line()

diff_info <- expand_grid(r=seq(0.02, 0.2, length.out=200), q0=seq(1e-5, 0.5, length.out=200)) |> 
    rowwise() |> 
    mutate(
        dom=dom_time(0, q0, r),
        crit=info_root_t(r, q0)$root
    )

## OK. A problem here is you've been thinking of the domination time as the symmetry point, 
## but this is NOT the case. They are very different as q0 -> 0.5 for example!!!

## ^^ 2/28 NOT true. You've shown the critical point for the derivative of the growth curve
## itself is indeed the domination time (proof: diff eq form of the model is rq(1 - q) which is
## always maximum at q=0.5)!

## thinking through below plot, it makes more sense now that for q0 -> 0.5, the inflection point
## is like t=1 but for low r the root will happen increasingly later

ggplot(diff_info, aes(r, q0)) +
    geom_raster(aes(fill=crit-dom), interpolate=TRUE) +
    scale_y_continuous(expand=expansion()) +
    scale_x_continuous(expand=expansion()) +
    scale_fill_viridis_c(option="mako", rescaler=~scales::rescale(.x)^0.5)

filter(diff_info, r > 0, q0 < 0.1, i2 < i1)


1 - nopt_single_param(tmax+1, 0, q0, r, error=0.15*r) / nopt_single_param(tmax, 0, q0, r, error=0.15*r)

    
# plots for desired prevalence on log-log scale-----------------------------------
q0 <- 1e-4
# qmax <- c(0.02, 0.1)
# r <- seq(0.05, 0.3, length.out=60)
qmax <- seq(0.01, 0.6, length.out=60)
r <- round(seq(0.02, 0.3, length.out=6), 2)

vgrid_qmax <- expand_grid(qmax, r) |>
    rowwise() |>
    mutate(
        tmax=find_tmax(qmax, 0, q0, r),
        nopt=nopt_single_param(tmax, 0, q0, r, error=r*0.15),
        nopt_tot=nopt*tmax, 
        var=var_single_param(nopt, tmax, 0, q0, r)
    )

lm(log(vgrid_qmax$nopt) ~ log(vgrid_qmax$qmax))

# percentage decrease in samples needed for a 1% increase in prevalence
(1.01^(-1.314) - 1) * 100

vgrid_qmax |>
    group_by(tmax, r) |>
    slice_min(qmax, n=1) |>
    ungroup() |>
    # mutate(label=paste0("q[T]=", qmax)) |>
    filter(nopt >= 1) |> 
    ggplot(aes(qmax, nopt, group=factor(r), col=factor(r))) +
    geom_line(linewidth=1.2) +
    # geom_abline(slope=-1.314, intercept=0) +
    scale_color_manual(values=pal_discrete[1:6]) +
    scale_y_continuous(transform="log", breaks=c(50, 100, 150, 200)) +
    scale_x_continuous(transform="log", breaks=c(0.01, seq(0.1, 0.5, 0.2)), labels=scales::label_percent()) +
    # scale_x_continuous(transform="logit", breaks=c(0.01, seq(0.1, 0.9, 0.2), 0.99)) +
    labs(x="prevalence threshold", y="sequences per day", col="growth\nrate") +
    theme_half_open() +
    theme(
        legend.position="inside",
        legend.position.inside=c(0.74, 0.74)
    )

p2 <- vgrid_qmax |>
    group_by(tmax, r) |>
    slice_min(qmax, n=1) |>
    ungroup() |>
    ggplot(aes(qmax, nopt_tot, group=factor(r), col=factor(r))) +
    geom_line(linewidth=1.2) +
    scale_color_manual(values=pal_discrete[1:6]) +
    scale_y_continuous(transform="log", breaks=c(100, 500, 1:4*1000)) +
    scale_x_continuous(transform="log", breaks=c(0.01, seq(0.1, 0.5, 0.2)), labels=scales::label_percent()) +
    # scale_x_continuous(transform="logit", breaks=c(0.01, seq(0.1, 0.9, 0.2), 0.99)) +
    labs(x="prevalence threshold", y="total sequences", col="growth\nrate") +
    theme_half_open() +
    theme(
        legend.position="inside",
        legend.position.inside=c(0.74, 0.74)
    )

plot_grid(p1, p2, nrow=1, rel_widths=c(1, 1.03))
ggsave("figs/nopt-single-param-qmax.pdf", width=7.5, height=3.75)

# key points:
# waiting until VOC has grown will significantly reduce the number of samples needed (per day)
# however, the total number of samples you will need for a particular prevalence threshold is pretty 
# much the same regardless of how fast the variant grows

#  variant proportion over time for single param model----------------------------
d_prev_single_param <- function(t, q0, r) {
    a <- 1/q0 - 1
    eterm <- exp(-r * t)
    denom <- (1 + a*eterm)^2
    dqdr <- a * t * eterm / denom
    return(dqdr)
}

d_dom_single_param <- function(t, q0, r) { # TODO unused t somewhat awkward
    a <- 1/q0 - 1
    -log(a) / r^2
}

d_tmax_single_param <- function(q, t, q0, r) {
    a <- 1/q0 - 1
    -log(-(q - 1)/(a*q))/r^2
}

nopt_g_single_param <- function(t, q0, r, dgdr, error=0.05, level=0.95, tmax=t) { # generally would assume we only observe up to t as well
    vr <- var_single_param(tmax, q0, r)
    vg <- vr * dgdr(t, q0, r)^2
    Z <- -qnorm((1 - level)/2)
    return((Z/error)^2 * vg)
}

nopt_prev_classic <- function(t, q0, r, error=0.05, level=0.95) {
    Z <- -qnorm((1 - level)/2)
    qtrue <- final_prevalence(t, 0, q0, r)
    (Z/error)^2 * qtrue * (1 - qtrue)
}



# for presentations: show rhat proogating uncertainty in prevalence
r <- 0.07 * 7 # weekly sequencing for this example
q0 <- 1e-2
qmax <- 0.48 # setting 0.5 will cause overshoot to 0.6 with these weekly settings
tmax <- find_tmax(qmax, 0, q0, r)
n <- 10

qtrue <- growth_model2(tmax, 0, q0, r)

vr <- var_single_param(tmax, q0, r) / n
vprev <- vr * d_prev_single_param(0:tmax, q0, r)^2

ggplot() +
    geom_area(aes(x=r+c(-0.25, 0.25)), fun=\(x) dnorm(x, r, sqrt(vr)), col="gray70", fill="#549DA7", stat="function") +
    # xlim(r+c(-0.15, 0.15)) +
    labs(x="growth rate", y="density") +
    theme_half_open()
    # theme(
    #     axis.text=element_blank(),
    #     axis.ticks=element_blank()
    # )

ggsave("figs/delta-method-ex-a.pdf", width=3.8, height=3.4)

tibble(t=0:tmax, q=qtrue, vprev) |> 
    ggplot(aes(t)) +
    geom_ribbon(aes(ymin=qnorm(0.025, q, sqrt(vprev)), ymax=qnorm(0.975, q, sqrt(vprev))), fill="#549DA7", col=NA) +
    geom_line(aes(y=q), col="black") +
    scale_x_continuous(breaks=c(0, 5, 10)) +
    labs(x="weeks", y="prevalence") +
    theme_half_open()

ggsave("figs/delta-method-ex-b.pdf", width=3.8, height=3.4)

# compare with classic version as function of sample size
r <- 0.07
q0 <- 1e-4
qmax <- seq(0.01, 0.95, length.out=100)

qmax_classic_comp <- tibble(tmax=0:180) |> 
    rowwise() |> 
    mutate(
        qmax=final_prevalence(tmax, 0, q0, r),
        `Classic result`=nopt_prev_classic(tmax, q0, r, error=qmax*0.15),
        `Delta method`=nopt_g_single_param(tmax, q0, r, d_prev_single_param, error=qmax*0.15)
    ) |> 
    ungroup()

p1 <- tibble(tmax=0:180, q=growth_model2(180, 0, q0, r)) |> 
    # slice(-1) |> 
    ggplot(aes(tmax, q)) +
    geom_line(col="gray50", linewidth=1.06) +
    geom_vline(xintercept=66, linetype="dashed", linewidth=1.04, col="gray50") +
    labs(x="days of observation", y="prevalence") +
    theme_half_open()

p2 <- qmax_classic_comp |> 
    filter(qmax > 0.01) |> 
    pivot_longer(c(`Classic result`, `Delta method`)) |> 
    ggplot(aes(tmax, value, col=name)) +
    geom_line(linewidth=1.06) +
    geom_vline(xintercept=66, linetype="dashed", linewidth=1.04, col="gray50") +
    scale_y_continuous(
        transform="log", breaks=10^(0:5)
        # sec.axis=sec_axis(
        #     transform=~./ymax,
        #     breaks=seq(0, 1, 0.2)
        # )
    ) +
    # scale_x_continuous(sec.axis=dup_axis(
    #     name="prevalence", 
    #     breaks=qmax_labs$b, 
    #     labels=scales::scientific(qmax_labs$l, 1),
    #     guide=guide_axis(angle=45)
    # )) +
    scale_color_manual(values=c(pal_discrete[c(1, 6)])) +
    xlim(c(0, 180)) +
    labs(x="days of observation", y="sequences per day", col=NULL) +
    theme_half_open() +
    theme(legend.position=c(0.52, 0.92))

plot_grid(p1, p2, nrow=1)
ggsave("figs/nopt-prev-classic-comp.pdf", width=7.7, height=3.5)

# sample size as func of prevalence for several metrics
# TODO: the domination/concern times just seem really off. Why would sample sizes be so big??
nopt_prev_qmax <- tibble(qmax) |>
    rowwise() |>
    mutate(
        tmax=find_tmax(qmax, 0, q0, r),
        dom_time=find_tmax(0.5, 0, q0, r),
        nopt_r=nopt_single_param(tmax, q0, r, error=r*0.15),
        nopt_prev=nopt_g_single_param(tmax, q0, r, d_prev_single_param, error=qmax*0.15),
        nopt_dom=nopt_g_single_param(tmax, q0, r, d_dom_single_param, error=1),
        nopt_concern=nopt_g_single_param(
            tmax, q0, r, 
            \(t, q0, r) d_tmax_single_param(0.05, t, q0, r), 
            error=1
        )
    )

nopt_prev_qmax |> 
    group_by(tmax) |>
    slice_min(qmax, n=1) |>
    ungroup() |>
    pivot_longer(contains("nopt")) |> 
    filter(value > 1/7, name %in% c("nopt_r", "nopt_prev")) |> 
    # filter(value > 1/7) |> 
    mutate(
        total=value*tmax,
        label=factor(fct_inorder(name), labels=c("growth rate", "prevalence"))
        # label=factor(fct_inorder(name), labels=c("growth rate", "prevalence", "domination time", "concern time"))
    ) |> 
    ggplot(aes(qmax, value, col=label)) +
    geom_line(linewidth=1.2) +
    scale_color_manual(values=pal_discrete[c(1, 3, 5, 7)]) +
    scale_y_continuous(transform="log10") +
    scale_x_continuous(transform="log10", breaks=c(0.01, 0.1, 0.5, 0.99), labels=scales::label_percent()) +
    # scale_x_continuous(transform="logit", breaks=c(0.01, seq(0.1, 0.9, 0.2), 0.99)) +
    labs(x="prevalence threshold", y="sequences per day", col=NULL) +
    theme_half_open()

ggsave("figs/nopt-growth-prev-compare.pdf", width=5.3, height=3.75)

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
