library(tidyverse)
library(lubridate)
library(jsonlite)
library(cowplot)

source("src/growth-model.R")
source("src/nopt-single-param.R")

# at date of first detection of Omicron (21L), 11494 sequences reported in Mass
# ^^^also the maximum sequencing effor I believe
# following weeks were 

voc0 <- fromJSON("data/perCountryData.json")

covid0 <- read_csv("data/nhsn-hospital-admissions.csv")

voc0$regions

voc_mass <- voc0$regions |> 
    as_tibble() |> 
    filter(region == "United States") |> 
    # unnest(cluster_names) |> 
    unnest(distributions) |> 
    filter(country == "Massachusetts") |> 
    unnest(distribution)

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
        nopt_easy=nopt_single_param(tmax, q0, r, error=r*0.25, level=0.75),
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
