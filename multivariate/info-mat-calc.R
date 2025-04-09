library(tidyverse)
library(cowplot)
library(MASS, exclude=c("select"))

source("src/growth-model.R")

# TODO: you can use an analytical gradient for the Monte Carlo method!

voc_growth_error <- function(n, cov, level=0.95) {
    Z <- -qnorm((1 - level)/2)
    Z * sqrt(cov[1, 1]/n)
}

# dumb algorithm that has the benefit that n>1 always errs on the too larger side, 
# within an error tolerance of `etol`
find_sample_size <- function(tmax, t0, q0, r, known_idx=NULL, level=0.95, error=0.05) {
    n <- 1
    # print(tmax)
    im <- round(info_mat(tmax, t0, q0, r), 6)
    idx <- 1:3
    if (!is.null(known_idx)) {
        idx <- idx[-known_idx]
    }
    cov <- solve(im[idx, idx])
    err <- voc_growth_error(n, cov, level)
    # print(cov)
    while (err > error) {
        n <- 2*n
        err <- voc_growth_error(n, cov, level)
    }
    while (err <= error) { # now walk backwards to an exact n
        n <- n-1
        err <- voc_growth_error(n, cov, level)
    }
    return(n+1)
}

# optimal sample size as function of timespan (tmax - t0 [assuming t0_true > 0])--
r_alpha <- 0.07 * 7
r_omicron <- 0.12 * 7
t0 <- 2.01
tmax <- 4:24

nopt_alpha <- map_dbl(tmax, ~find_sample_size(., t0, 1e-4, r_alpha))

nopt_res_alpha <- tibble(
    tmax, 
    nopt_week=nopt_alpha, nopt_total=tmax * nopt_alpha,
    qmax=map_dbl(tmax, ~final_prevalence(., t0, 1e-4, r_alpha))
)

nopt_res_omicron <- tibble(
    tmax, 
    nopt_week=nopt_omicron, nopt_total=tmax * nopt_omicron,
    qmax=map_dbl(tmax, ~final_prevalence(., t0, 1e-4, r_omicron))
)

nopt_res <- bind_rows(
    mutate(nopt_res_alpha, variant="alpha"),
    mutate(nopt_res_omicron, variant="omicron")
)

scale_nopt <- function(n, l, u) {
    (log10(n) - l)/(u-l)
}

inv_scale_nopt <- function(x, l, u) {
    10^(x*(u-l) + l)
}

upper <- log10(1e8)
lower <- 0.95*log10(min(nopt_res$nopt_total))

ggplot(nopt_res, aes(tmax-2)) +
    # geom_line(aes(y=tmax*max(tot_seq$n)/2), col="#7fcdbb", linewidth=1.1) +
    # geom_vline(xintercept=6.7, col="#7fcdbb", linewidth=1.05) +
    # geom_line(aes(y=tmax*median(tot_seq$n)/2), col="#1d91c0", linewidth=1.1) +
    # geom_vline(xintercept=9.95, col="#1d91c0", linewidth=1.05) +
    geom_hline(
        yintercept=scale_nopt(10^(4:8), lower, upper), 
        linetype=32, col="white", linewidth=1
    ) +
    geom_line(aes(y=qmax, col=variant), linewidth=1.1, linetype=112) +
    geom_line(aes(y=scale_nopt(nopt_total, lower, upper), col=variant), linewidth=1.13) +
    geom_point(aes(y=scale_nopt(nopt_total, lower, upper), col=variant), size=1.4) +
    # facet_wrap(~name, scales="free", nrow=2) +
    scale_color_manual(values=c("darkorange1", "#35A0ABFF")) +
    scale_y_continuous(
        position="right", limits=c(0, 1),
        sec.axis=sec_axis(
            ~inv_scale_nopt(., lower, upper), 
            name="total sequences needed", breaks=10^(4:8), labels=scales::label_comma()(10^(4:8))
        )
    ) +
    # scale_y_log10(
    #     labels=scales::comma
    #     # sec.axis=sec_axis(~./sum(.))
    # ) +
    # scale_x_continuous(
    #     breaks=tmax-2,
    #     sec.axis=dup_axis(
    #         name="final prevalence", 
    #         labels=str_remove(scales::scientific(nopt_res_alpha$qmax, 1, trim=FALSE), "0"),
    #         guide=guide_axis(angle=45)
    #     ),
    # ) +
    labs(x="weeks since introduction", y="final prevalence", linetype=NULL) +
    theme_minimal_grid() +
    background_grid(major="none") +
    theme(
        legend.position="none", 
        axis.ticks=element_line(linewidth=1),
        panel.background=element_rect(fill="gray90")
    )

ggsave("figs/nopt-omicron.pdf", width=5.5, height=3.2)

# minimum total sequences needed for different growth rates-----------------------
find_min_tmax <- function(t0, r, q0=1e-4) {
    tmax <- 5
    nopt <- nopt_prev <- tmax*find_sample_size(tmax, t0, q0, r)
    while (nopt <= nopt_prev) {
        nopt_prev <- nopt
        tmax <- tmax+1
        nopt <- tmax*find_sample_size(tmax, t0, q0, r)
    }
    return(tibble_row(r=r, tmax=tmax-1, n=nopt_prev))
}

t0 <- 2.01
r <- seq(0.05, 1.5, length.out=25)

tmax_res <- map_dfr(r, ~find_min_tmax(t0, .))

tmax_res |> 
    distinct(tmax, .keep_all=TRUE) |> 
    ggplot(aes(tmax, n, col=r)) +
    geom_line(linewidth=1.25) +
    scale_color_viridis_c() +
    theme_minimal_grid() +
    background_grid(minor="none") +
    theme(
        axis.ticks=element_line(linewidth=1),
        panel.background=element_rect(fill="gray90"),
        panel.grid=element_line(color="white", linewidth=1)
    )

# final prevalence using delta method---------------------------------------------

# TODO: this classic version throws out all data before tmax.
# question is, can you get smaller CIs sometimes by using all data and calculating
# q_mle = f(mle) ?
voc_fin_prev_err_classic <- function(n, tmax, t0, q0=1e-4, r, level=0.95) {
    qmax <- final_prevalence(tmax, t0, q0, r)
    var <- (qmax*(1-qmax)) / n
    1 / im
    ....
}
# TODO its definitely time soon to cleanup this script's organization and naming
# - remove voc_* naming
# - reduntant find sample size functions
# - split into methods/helpers and result plotting

grad_fin_prev <- function(tmax, t0, r, q0) {
    a <- 1/q0 - 1
    eterm <- exp(-r * (tmax-t0))
    denom <- (1 + a*eterm)^2
    c(a * (tmax-t0) * eterm / denom, -a * r * eterm / denom)
}

voc_fin_prev_error <- function(n, cov, grad_tmax, level=0.95) {
    var_qmax <- t(grad_tmax) %*% (cov/n) %*% grad_tmax
    Z <- -qnorm((1 - level)/2)
    Z * sqrt(var_qmax)
}

find_nopt_fin_prev <- function(tmax, t0, q0, r, level=0.95, error=0.05) {
    # n <- 1
    # cov <- solve(info_mat(tmax, t0, q0, r))
    # err <- voc_growth_error(n, cov, level)
    # while (err > error) {
    #     n <- 2*n
    #     err <- voc_growth_error(n, cov, level)
    # }
    # while (err <= error) { # now walk backwards to an exact n
    #     n <- n-1
    #     err <- voc_growth_error(n, cov, level)
    # }
    # return(n+1)
}

# just was for confirming the ugly equation version matched matrix version
voc_fin_prev_err3 <- function(n, tmax, t0, q0=1e-4, r, level=0.95) {
    im <- info_mat(n, tmax, t0, q0=1e-4, r)
    cv <- solve(im)
    a <- 1/q0 - 1
    eterm <- exp(-r * (tmax-t0))
    fT <- 1 + a*eterm
    dfTdr <- a*(tmax-t0)*eterm
    dfTdt0 <- -a*r*eterm
    (dfTdr^2*cv[1, 1] + dfTdt0^2*cv[2, 2] + 2*dfTdr*dfTdt0*cv[1, 2]) / fT^4
}

r_alpha <- 0.07 * 7
t0 <- 2.01
tmax <- 30

cov <- solve(info_mat(tmax, t0, 1e-4, r_alpha))
jqmax <- grad_fin_prev(tmax, t0, r_alpha, 1e-4)
voc_fin_prev_error(1000, cov, jqmax)
voc_growth_error(1000, cov)

cov <- solve(info_mat(1000, tmax, t0, 1e-4, r_alpha))
mle_approx <- MASS::mvrnorm(100000, c(r_alpha, t0), cov)

-qnorm(0.05/2) * sqrt(var(map2_dbl(mle_approx[,1], mle_approx[,2], \(r, t0) final_prevalence(tmax, t0, 1e-4, r))))
