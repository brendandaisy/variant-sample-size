library(tidyverse)
library(lubridate)
library(rjson)

voc0 <- fromJSON(file="data/perCountryData.json")
length(voc0$regions[[2]]$distributions)
voc_us <- voc$regions[[1]]$distributions[[1]]$distribution # first [[1]] is to get "World" tab, second is to get all "USA" data

voc_us_df <- map_dfr(voc_us, \(x) {
    counts <- as_vector(x$cluster_counts)
    counts <- c(counts, others=(x$total_sequences - sum(counts)))
    ret <- enframe(counts, value="count")
    
    mutate(ret, total_sequences=x$total_sequences, prop=counts/total_sequences, week=ymd(x$week))
})

voc_us_df |> 
    group_by(week, total_sequences) |> 
    reframe(count=sum(count))

p2 <- voc_us_df |> 
    filter(week > "2022-11-01", week < "2024-09-01") |> 
    distinct(total_sequences, week) |> 
    ggplot(aes(week, total_sequences)) +
    geom_col() +
    labs(y="sequences submitted") +
    theme_half_open()

p1 <- voc_us_df |> 
    filter(week > "2022-11-01", week < "2024-09-01", count > 10) |> 
    ggplot(aes(week, prop)) +
    geom_area(aes(fill=name), position=position_fill()) +
    labs(x="week", y="observed proportion", fill="nextclade\nvariant") +
    scale_x_date(expand=expansion()) +
    scale_y_continuous(expand=expansion()) +
    theme_half_open()

plot_grid(p1, p2, nrow=1, rel_widths=c(1, 0.7))
ggsave("figs/covariants-data.pdf", width=10.2, height=4.8)


### compare nopt for Omicron like variant and state sequencing activity
r <- 0.12
q0 <- 1e-4
t0 <- 0

qmax <- seq(0.01, 0.75, length.out=60)

nopt_omicron <- tibble(qmax) |>
    rowwise() |>
    mutate(
        tmax=find_tmax(qmax, t0, q0, r),
        nopt=nopt_single_param(tmax, t0, q0, r, error=r*0.15),
        nopt_tot=nopt*tmax, 
    ) |> 
    group_by(tmax) |>
    slice_min(qmax, n=1) |>
    ungroup()

n_state <- voc_us_df |> 
    group_by(year=year(week)) |> 
    summarise(best_eff=quantile(total_sequences / (50*14), 0.75))

n_state_tot <- nopt_omicron |> 
    expand_grid(n_state) |> 
    # mutate(neff=tmax*best_eff) |> 
    distinct(qmax, year, best_eff)

nopt_omicron |> 
    ggplot(aes(qmax, nopt)) +
    geom_line(aes(qmax, best_eff, color=as.factor(year)), n_state_tot, linewidth=1.2, linetype="dashed") +
    geom_line(linewidth=1.2, col="gray50") +
    scale_color_manual(values=pal_discrete[1:6]) +
    scale_y_continuous(transform="log", breaks=c(1, 10, 100)) +
    scale_x_continuous(transform="log", breaks=c(0.01, 0.1, 0.3, 0.7), labels=scales::label_percent()) +
    labs(x="prevalence threshold", y="sequences per day", col=NULL) +
    theme_half_open()

ggsave("figs/nopt-omicron-state-compare.pdf", width=4.8, height=3.75)
