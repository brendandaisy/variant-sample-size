library(tidyverse)
library(lubridate)
# library(jsonlite)
library(cowplot)
library(usmap)
library(sf)
library(ggforce)

source("src/growth-model.R")
source("src/nopt-single-param.R")

ma0 <- read_csv("data/tabula-7AUG25_Total sequences_bycounty_forUGA.csv")
ma_counties <- usmap::us_map(regions=c("counties"), include=c("MA")) |> 
    st_transform("WGS84")

small_counties <- c("Dukes County", "Barnstable County", "Nantucket County")

cities <- st_sf(
    name=c("Boston", "Springfield"),
    geom=st_sfc(
        st_point(c(-71.057778, 42.360278)),
        st_point(c(-72.590278, 42.101389)),
        crs="WGS84"
    ))

ma_counties <- tibble(
    county="Dukes/Nantucket/Barnstable counties",
    geom=st_union(filter(ma_counties, county %in% small_counties))
) |> 
    st_as_sf() |> 
    bind_rows(filter(ma_counties, !county %in% small_counties))

ma <- ma0 |> 
    mutate(
        year=parse_integer(str_extract(`MMWR Week`, "\\d{4}")),
        epiweek=parse_integer(str_extract(`MMWR Week`, "\\d{2}$"))
    ) |> 
    select(year, epiweek, county=County, seq_type=`Sequence Type`, n=`Number Sequenced`)

ma_seq_success <- ma |> 
    group_by(county, seq_type) |> 
    summarize(n=sum(n), .groups="drop") |> 
    pivot_wider(values_from=n, names_from=seq_type) |> 
    transmute(
        county, 
        prop_success=`Successfully sequenced`/(`Invalid/Unsatisfactory`+`Successfully sequenced`)
    )

ma_seq_success
summarize(ma_seq_success, mean(prop_success), sd(prop_success))

ma_avg_dec21 <- ma |> 
    filter(seq_type == "Successfully sequenced") |> 
    filter(year == 2021, epiweek >= epiweek("2021-12-01")) |> 
    group_by(county) |> 
    summarize(n_daily=mean(n, na.rm=TRUE)/7)

ma_avg_25 <- ma |> 
    filter(seq_type == "Successfully sequenced") |> 
    filter(year == 2025, epiweek >= epiweek("2025-01-01"), epiweek < epiweek("2025-07-01")) |> 
    group_by(county) |> 
    summarize(n_daily=mean(n, na.rm=TRUE)/7)

# State-level sequencing activity-------------------------------------------------
ma_state <- ma |> 
    filter(seq_type == "Successfully sequenced") |> 
    group_by(year, epiweek) |> 
    reframe(n=sum(n))

covid0 <- read_csv("data/nhsn-hospital-admissions-oct-8.csv")
first_wk <- min(filter(ma_state, year == 2021)$epiweek)
last_wk <- max(filter(ma_state, year == 2025)$epiweek)

mass_q0 <- covid0 |> 
    select(
        date=`Week Ending Date`,
        location=`Geographic aggregation`,
        cases=`Total Patients Hospitalized with COVID-19`
    ) |> 
    filter(
        location == "MA", 
        year(date) < 2025 | epiweek(date) < last_wk, 
        year(date) > 2021 | (year(date) == 2021 & epiweek(date) > first_wk & epiweek(date) != 53)
    ) |> 
    transmute(year=year(date), epiweek=epiweek(date), date, cases, q0=1/cases)

# start with a daily growth rate, compute nopt at the start of each week, then 
# convert back to weekly samples
r <- 0.84
tmax <- 4 # let's just use tmax since qmax leads to jumps when gaining a week

nopt_state_r <- mass_q0 |> 
    filter(is.finite(q0), q0 < 0.05) |> # weeks with <20 reported cases removed
    rowwise() |> 
    mutate(
        # tmax_daily=max(find_tmax(qmax, q0, r), tmax_min),
        tmax_daily=tmax*7,
        nopt_easy=nopt_single_param(tmax, q0, r, error=r*0.25, level=0.2),
        nopt_hard=nopt_single_param(tmax, q0, r, error=r*0.15)
    ) |> 
    ungroup() |> 
    complete(date=unique(mass_q0$date)) |> 
    pivot_longer(contains("nopt"))

ggplot(mass_q0, aes(date, cases)) +
    geom_line()

ma_state |> 
    left_join(mass_q0) |> 
    ggplot(aes(date, n)) + #convert to weekly scale
    geom_col(aes(fill="A"), col="white", linewidth=0.05) +
    geom_line(aes(date, value, col=name), nopt_state_r, linewidth=1.05) +
    scale_y_continuous(
        transform="log10", breaks=10^(1:4)
        # sec.axis=dup_axis(name="required sequences per week", labels=NULL)
    ) +
    scale_fill_manual(values=c("gray50"), labels=c("reported sequences")) +
    scale_color_manual(
        values=pal_discrete[c(1, 6)],
        labels=c("required sequences (low accuracy)", "required sequences (high accuracy)")
    ) +
    # scale_x_date(date_breaks="6 months", date_labels="%m-%y", guide=guide_axis(angle=45)) +
    labs(x="week", y="sample size", col=NULL, fill=NULL) +
    theme_half_open() +
    theme(
        # legend.position="inside",
        # legend.position.inside=c(0.4, 0.85), 
        # axis.ticks.y.right=element_blank(),
        legend.text=element_text(size=rel(0.7)),
        legend.spacing=unit(0.5, "mm")
    )

ggsave("figs/mass-state-level-2.pdf", width=7.4, height=4)


# TODO: don't really love this. Let's try prevalence in the next figure and do something
# more positive / high level with this panel
nopt_state_prev <- mass_q0 |> 
    filter(is.finite(q0), q0 < 0.05) |> # weeks with <20 reported cases removed
    rowwise() |> 
    mutate(
        tmax_daily=max(find_tmax(qmax, q0, r), tmax_min),
        nopt_easy_prev=nopt_g_single_param(
            tmax_daily, q0, r, d_prev_single_param, error=qmax*0.3, level=0.2
        )*7,
        nopt_hard_prev=nopt_g_single_param(tmax_daily, q0, r, d_prev_single_param, error=qmax*0.15)*7
    ) |> 
    ungroup() |> 
    complete(date=unique(mass_q0$date)) |> 
    pivot_longer(contains("nopt"))

slice_max(nopt_state_prev, value)

ma_state |> 
    left_join(mass_q0) |> 
    ggplot(aes(date, n)) + #convert to weekly scale
    geom_col(aes(fill="A"), col="white", linewidth=0.1) +
    geom_line(aes(date, value, col=name), nopt_state_prev, linewidth=1.05) +
    scale_y_continuous(transform="log10", breaks=10^(1:4)) +
    scale_fill_manual(values=c("gray50"), labels=c("reported sequences")) +
    scale_color_manual(
        values=pal_discrete[c(1, 6)],
        labels=c("required sequences (low accuracy)", "required sequences (high accuracy)")
    ) +
    # scale_x_date(date_breaks="6 months", date_labels="%m-%y", guide=guide_axis(angle=45)) +
    labs(x="week", y="reported sequences", col=NULL, fill=NULL) +
    theme_half_open() +
    theme(legend.position="none", axis.ticks.y.right=element_blank())

# How long after it reached 1% prevalence would it have taken a county----------
# to detect an Omicron-like variant at different points?--------------------------
# How about to be able to track prevalence?---------------------------------------
q0 <- 1e-2
r <- 0.12

ma_tmin_dec21 <- ma_avg_dec21 |>
    rowwise() |> 
    mutate(
        tmin_r=nopt_tmin(n_daily, q0, r, level=0.05),
        tmin_prev=nopt_tmin(n_daily, q0, r, do_prev=TRUE, error=0.15, level=0.05)
    ) |> 
    ungroup()

ma_tmin_25 <- ma_avg_25 |>
    rowwise() |> 
    mutate(
        tmin_r=nopt_tmin(n_daily, q0, r, level=0.05),
        tmin_prev=nopt_tmin(n_daily, q0, r, do_prev=TRUE, error=0.15, level=0.05)
    ) |> 
    ungroup()

ma_tmin <- bind_rows(
    mutate(ma_tmin_dec21, period="December 2021"), 
    mutate(ma_tmin_25, period="January-June 2025")
)

curve <- tibble(
    t=0:65,
    q=growth_model(65, q0, r)
)

p1 <- curve |> 
    ggplot(aes(t, q)) +
    geom_line(
        data=filter(curve, t > max(ma_tmin$tmin_r, na.rm=TRUE)),
        col="gray50", linewidth=2, linetype="11"
    ) +
    geom_line(
        data=filter(curve, t <= min(ma_tmin$tmin_r, na.rm=TRUE)),
        col="black", linewidth=2, linetype="11"
    ) +
    geom_line(
        aes(col=t), 
        filter(curve, t >= min(ma_tmin$tmin_r, na.rm=TRUE), t <= max(ma_tmin$tmin_r, na.rm=TRUE)), 
        linewidth=3.4, lineend="square"
    ) +
    annotate(
        "label", 61, 0.6,
        label="unable to\ndetect", col="gray50", label.size=NA
    ) +
    scale_color_gradientn(colors=pal_discrete) +
    scale_y_continuous(labels=scales::label_percent()) +
    # scale_x_continuous(sec.axis=dup_axis(name="days needed")) +
    labs(x="days needed", y="prevalence") +
    theme_half_open() +
    theme(
        legend.position="none", 
        panel.grid.major.y=element_line(linetype="dashed", color="gray70"),
        plot.margin=margin(0.1, 0.82, 0, 0.65, "in")
        # axis.text.x.bottom=element_blank(),
        # axis.title.x.bottom=element_blank(),
        # axis.line.x.bottom=element_blank(),
        # axis.ticks.x.bottom=element_blank()
    )

p2 <- ma_counties |> 
    left_join(ma_tmin) |> 
    pivot_longer(contains("tmin"), names_to="target", values_to="tmin") |> 
    mutate(target=fct_inorder(factor(target, labels=c("prevalence", "growth advantage")))) |> 
    ggplot() +
    geom_sf(aes(fill=tmin), col="white") +
    geom_sf(data=cities) +
    geom_text(
        aes(X, Y, label=name), bind_cols(st_drop_geometry(cities), st_coordinates(cities)),
        nudge_x=c(0.43, -0.2), nudge_y=c(0, -0.12)
    ) +
    facet_grid(target~period) +
    scale_fill_gradientn(colors=pal_discrete, na.value="gray50") +
    labs(fill="days needed") +
    theme_map() +
    theme(
        strip.text=element_text(face="bold"),
        legend.position="none",
        plot.margin=margin(0, 0, 0, 0),
        # legend.position="bottom",
        # legend.justification="center",
        # legend.key.width=unit(1, "in"),
        # legend.title.position="top",
        # legend.title=element_text(hjust=0.5)
    )

# ggplot(tibble(x=unique(ma_tmin$tmin)), aes(x, y=0, col=x)) +
#     geom_line(linewidth=3) +
#     scale_color_gradientn(colors=pal_discrete) +
#     scale_x_continuous(sec.axis=dup_axis())
#     labs(x="days needed") +
#     theme_half_open() +
#     theme(
#         legend.position="none",
#         axis.line.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         axis.title.y=element_blank()
#     )

# plot_grid(plot_grid(
    # p1, plot_grid(NULL, p2, NULL, rel_widths=c(0.1, 0.8, 0.1), nrow=1), 
#     p1, p2,
#     rel_heights=c(1, 1),
#     ncol=1
# ),
# legend, nrow=1, rel_widths=c(1, 0.12))

plot_grid(p1, p2, rel_heights=c(0.75, 2), ncol=1, labels="AUTO")

ggsave("figs/ma-county-tmin-3.pdf", width=5.6, height=5.8)

# Sample size to weight hypotheses of lower growth following introduction---------
# in other regions----------------------------------------------------------------

# ma_t: ma sequences counts for a single time unit (1 row per county)
nopt_growth_compare_sit1 <- function(ma_t, q0A, q0, rA, tmax, pct_dec=0.15, alpha=0.05, beta=0.2,
                                     countyA="Suffolk County", just_tstat=FALSE) {
    M <- nrow(ma_t)
    n_prop <- ma_t$n_daily / sum(ma_t$n_daily)
    
    r <- rep(rA - rA*pct_dec, M)
    r[which(ma_t$county == countyA)] <- rA
    q0 <- rep(q0, M)
    q0[which(ma_t$county == countyA)] <- q0A
    
    if (just_tstat) {
        return(wald_stat(tmax, q0, r, M))
    }
    
    nopt_multiple_regions(tmax, q0, r, M, n_prop, 1-beta, alpha)
}

power_growth_compare_sit1 <- function(ma_t, q0A, q0, rA, tmax, pct_dec=0.15, countyA="Suffolk County") {
    M <- nrow(ma_t)
    n <- ma_t$n

    if (M <= 1 | !(countyA %in% ma_t$county))
        return(NA)
    
    r <- rep(rA - rA*pct_dec, M)
    r[which(ma_t$county == countyA)] <- rA
    q0 <- rep(q0, M)
    q0[which(ma_t$county == countyA)] <- q0A
    
    power_multiple_regions(tmax, q0, r, M, n)
}

nopt_dec21_t1t2 <- expand_grid(alpha=c(0.001, 0.01, 0.05), beta=c(0.01, 0.05, 0.2)) |> 
    rowwise() |> 
    mutate(
        nopt=nopt_growth_compare_sit1(
            ma_avg_dec21, 0.05, 1e-3, 0.12, 30,
            alpha=alpha, beta=beta, countyA="Suffolk County"
        )
    ) |> 
    ungroup() |> 
    mutate(alpha=fct_inorder(str_c(alpha*100, "%")), beta=fct_inorder(str_c(beta*100, "%")))

ggplot(nopt_dec21_t1t2, aes(alpha, beta)) +
    geom_tile(aes(fill=nopt)) +
    geom_text(aes(label=floor(nopt)), col="gray30") +
    scale_fill_gradientn(colors=pal_discrete) +
    scale_x_discrete(expand=expansion()) +
    scale_y_discrete(expand=expansion()) +
    labs(x="type I error", y="type II error") +
    theme_half_open() +
    theme(
        legend.position="none",
        axis.line=element_blank(),
        text=element_text(color="gray30"),
        axis.text=element_text(color="gray30")
    )

nopt_dec21_geo <- tibble(county=unique(ma$county)) |> 
    rowwise() |> 
    mutate(
        nopt=nopt_growth_compare_sit1(
            ma_avg_dec21, 0.05, 1e-3, 0.12, 30,
            countyA=county
        )
    ) |> 
    ungroup()

ma_counties |> 
    left_join(nopt_dec21_geo) |> 
    ggplot() +
    geom_sf(aes(fill=nopt), col="white") +
    scale_fill_gradientn(colors=pal_discrete) +
    # labs(fill="days needed", title="Time until an Omicron-like growth advantage could be\ninferred with 15% accuracy, 80% of the time") +
    theme_map()

### Type 1 and 2 error probabilities by week under different proportions:

# 9/8: of course alpha is a good match simulating, its baked right into `level` !
# sim_err_multiple_regions <- function(tmax, q0, r_null, r_alt, M, n, level=0.05, nsim=1000) {
#     crit_t <- qchisq(1-level, M-1) # critical value we need tstat to exceed to reject null
#     
#     var_null <- map_dbl(q0, ~var_single_param(tmax, .x, r_null)) / n
#     var_alt <- map2_dbl(q0, r_alt, ~var_single_param(tmax, .x, .y)) / n
#     
#     alpha <- map_lgl(1:nsim, \(i) {
#         r0_samp <- map_dbl(var_null, ~rnorm(1, r_null, sqrt(.x)))
#         vc <- diag(map2_dbl(q0, r0_samp, ~var_single_param(tmax, .x, .y)) / n)
#         C <- cbind(rep(-1, M-1), diag(1, M-1))
#         W <- solve(C %*% vc %*% t(C))
#         tstat <- t(C %*% r0_samp) %*% W %*% (C %*% r0_samp)
#         tstat > crit_t # accidentally reject
#     }) |> mean()
#     
#     beta <- map_lgl(1:nsim, \(i) {
#         r1_samp <- map2_dbl(r_alt, var_alt, ~rnorm(1, .x, sqrt(.y)))
#         vc <- diag(map2_dbl(q0, r1_samp, ~var_single_param(tmax, .x, .y)) / n)
#         C <- cbind(rep(-1, M-1), diag(1, M-1))
#         W <- solve(C %*% vc %*% t(C))
#         tstat <- t(C %*% r1_samp) %*% W %*% (C %*% r1_samp)
#         tstat < crit_t # fail to reject
#     }) |> mean()
#     
#     list(alpha_sim=alpha, beta_sim=beta)
# }

rA <- 0.12*7
tmax <- 4
q0A <- 0.01
q0 <- 1e-3

ex_tstat50 <- nopt_growth_compare_sit1(ma_avg_dec21, q0A, q0, rA, tmax, pct_dec=0.50, just_tstat=TRUE)
ex_tstat95 <- nopt_growth_compare_sit1(ma_avg_dec21, q0A, q0, rA, tmax, pct_dec=0.95, just_tstat=TRUE)
ex_crit_t <- qchisq(0.95, 11)

ggplot() +
    stat_function(
        fun=~dchisq(.x, 11, ncp=100*ex_tstat50), geom="area",
        col=NA, fill=pal_discrete[2],  xlim=c(ex_crit_t, 60)
    ) +
    stat_function(
        fun=~dchisq(.x, 11), geom="area",
        col=NA, fill=pal_discrete[5], xlim=c(0, ex_crit_t)
    ) +
    stat_function(fun=~dchisq(.x, 11), col=pal_discrete[7], linewidth=1.05) +
    stat_function(fun=~dchisq(.x, 11, ncp=100*ex_tstat50), col=pal_discrete[1], linewidth=1.05) +
    
    labs(x="test statistic", y="density") +
    xlim(0, 60) +
    theme_half_open() +
    theme(
        axis.text=element_blank(),
        axis.ticks=element_blank()
    )

ggsave("figs/type-1-2-error-ex.pdf", width=3.4, height=2.9)

pow_sit1 <- ma |> 
    filter(seq_type == "Successfully sequenced") |> 
    group_by(year, epiweek) |> 
    group_modify(\(ma_t, g) {
        tibble_row(
            date=as.Date(str_c(g$year, "-", g$epiweek, "-0"), "%Y-%U-%w"),
            power50=power_growth_compare_sit1(ma_t, q0A, q0, rA, tmax, pct_dec=0.50),
            power95=power_growth_compare_sit1(ma_t, q0A, q0, rA, tmax, pct_dec=0.95),
            n_tot=sum(ma_t$n)
        )
    }) |> 
    ungroup()

## For text: increasing to 95% reduction barely increases power
pow_sit1 |> 
    transmute(pinc=power95-power50) |> 
    summarize(mean(pinc, na.rm=TRUE))

# max_n_tot <- max(log(pow_sit1$n_tot))

pmain <- ggplot(pow_sit1, aes(date)) +
    # geom_col(aes(y=log(n_tot) / max_n_tot), col=NA, fill="gray50") +
    geom_line(aes(y=power), col=pal_discrete[2], linewidth=1.1) +
    # scale_y_continuous(sec.axis=sec_axis(
    #     ~exp(.x*max_n_tot), "total sequences",
    #     breaks=c(1, 10, 100, 1000, 8000)
    # )) +
    labs(x="weeks", y="power") +
    theme_half_open()

pinset <- ggplot() +
    geom_sf(
        aes(fill=issuff), mutate(ma_counties, issuff=county == "Suffolk County"),
        col="white"
    ) +
    scale_fill_manual(values=c("gray50", "tomato1")) +
    theme_map() +
    theme(legend.position="none")

ggdraw(pmain) +
    draw_plot(pinset, 0.53, 0.52, 0.42, 0.38)

ggsave("figs/mult-regions-sit-1.pdf", width=4.8, height=3.8)

# voc0 <- fromJSON("data/perCountryData.json")
# voc0$regions
# 
# voc_mass <- voc0$regions |> 
#     as_tibble() |> 
#     filter(region == "United States") |> 
#     # unnest(cluster_names) |> 
#     unnest(distributions) |> 
#     filter(country == "Massachusetts") |> 
#     unnest(distribution) |> 
#     mutate(week=as.Date(week))
# 
# save(voc_mass, file="data/voc-mass.RData")

data("voc-mass")
covid0 <- read_csv("data/nhsn-hospital-admissions.csv")

voc_mass_long <- voc_mass |> 
    unnest(cluster_counts) |> 
    select(week, total_sequences, `20A.EU2`:`S:677P.Pelican`) |>
    pivot_longer(`20A.EU2`:`S:677P.Pelican`)

voc_mass_long |> 
    filter(str_detect(name, "Delta")) |> 
    ggplot(aes(week, value, fill=name, col=name)) +
    geom_col(position=position_stack())
    # theme(legend.position="none")

variant_periods <- tribble(
    ~name, ~wmin, ~wmax,
    "21J (Delta)", 0, 36, 
    "21K (Omicron)", 0, 10
)

voc_proc <- voc_mass_long |> 
    filter(name %in% variant_periods$name, value > 0) |> 
    group_by(name) |> 
    arrange(week) |> 
    mutate(weeks_passed=2*(1:n()-1)) |> 
    left_join(variant_periods) |> 
    filter(weeks_passed <= wmax)

ggplot(voc_proc, aes(weeks_passed, value/total_sequences, col=name)) +
    geom_point(aes(size=total_sequences))

log_likelihood <- function(r, ys, ns, tmax, q0) {
    q <- growth_model(tmax, q0, r)
    sum(map_dbl(1:(tmax+1), \(i) dbinom(ys[i], ns[i], q[i], log=TRUE)))
}

r_mle <- voc_proc |> 
    group_modify(\(sdf, key) {
        tmax <- sdf$wmax[1]/2
        q0 <- sdf$value[1]/sdf$total_sequences[1]
        fun <- function(r) log_likelihood(r, sdf$value, sdf$total_sequences, tmax, q0)
        tibble_row(
            var=key, q0=q0, tmax=sdf$wmax[1],
            r=optimize(fun, c(0, 3.5), maximum=TRUE)$maximum/2
        )
    })

r_mle |> 
    rowwise() |> 
    mutate(q=list(growth_model(tmax, q0, r)), t=list(0:tmax)) |> 
    unnest(c(q, t)) |> 
    ggplot(aes(t, q, col=as.factor(r))) +
    geom_line() +
    geom_point(
        aes(weeks_passed, value/total_sequences, col=name, size=total_sequences),
        voc_proc
    ) +
    scale_color_manual(values=pal_discrete[c(1, 6, 1, 6)])

mass_tot_seq <- select(voc_mass, week, total_sequences) |> 
    mutate(week=ymd(week))

mass_q0 <- covid0 |> 
    select(
        date=`Week Ending Date`,
        location=`Geographic aggregation`,
        cases=`Total Patients Hospitalized with COVID-19`
    ) |> 
    filter(location == "MA") |> 
    transmute(date, cases, q0=1/cases)

r <- 0.12*7

nopt_mass <- mass_q0 |> 
    # filter(q0 < 1, is.finite(q0)) |> 
    rowwise() |> 
    mutate(
        tmax=4,
        nopt_easy=nopt_single_param(tmax, q0, r, error=r*0.25, level=0.25),
        nopt_hard=nopt_single_param(tmax, q0, r, error=r*0.15)
    ) |> 
    ungroup() |> 
    pivot_longer(contains("nopt")) |> 
    mutate(value=ifelse(!is.finite(value), NA, value)) # replace date with 0 cases with NA

ggplot(mass_tot_seq, aes(week, total_sequences/2)) + #convert to weekly scale
    geom_col(fill="gray50", col="white", linewidth=0.1) +
    geom_line(aes(date, value, col=name), nopt_mass, linewidth=1.05) +
    scale_y_continuous(
        transform="log10", breaks=10^(1:4),
        sec.axis=dup_axis(name="required sequences per week", labels=NULL)
    ) +
    scale_color_manual(
        values=pal_discrete[c(1, 6)],
        labels=c("low accuracy", "high accuracy")
    ) +
    # scale_x_date(date_breaks="6 months", date_labels="%m-%y", guide=guide_axis(angle=45)) +
    labs(x="week", y="reported sequences", col=NULL) +
    theme_half_open() +
    theme(legend.position=c(0.6, 0.85), axis.ticks.y.right=element_blank())

ggsave("figs/mass-state-level.pdf", width=5.4, height=4.1)

as_tibble(voc0)

names(voc0$regions)
length(voc0$regions[[2]]$distributions)
voc_us <- voc$regions[[1]]$distributions[[1]]$distribution
