library(tidyverse)

source("src/growth-model.R")

N0 <- 500 # total infectious cases when surveillance begins (not nec. at t0)
RH <- 0.2 # growth rate (beta - alpha, if you please) of infectious population for the historical strain
RV <- 0.6 # growth of VOC cases

t0 <- 3
tmax <- 12

res <- tibble(
    t=seq(0, tmax, 0.1), 
    V=ifelse(t<t0, NA, exp(RV*(t-t0))), 
    H1=(15-1) * exp(0.14*t), 
    H2=(15-1) * exp(-0.03*t), 
    H3=(90-1) * exp(0.085*t), 
    H4=(90-1) * exp(-0.03*t), 
)

p1 <- res |> 
    pivot_longer(-t) |> 
    ggplot(aes(t, value, col=name)) +
    geom_line() +
    scale_color_manual(values=c("#f5b43d", "#ff856d", "#c63df5", "#f53dbc", "black")) +
    labs(title="VOC growth compared to current trends", col="variant", x="weeks", y="case count") +
    theme_half_open()

p2 <- res |> 
    pivot_longer(H1:H4) |> 
    # group_by(name) |> 
    mutate(q=V/(V+value)) |> 
    ggplot(aes(t, q, col=name)) +
    geom_line() +
    scale_color_manual(values=c("#f5b43d", "#ff856d", "#c63df5", "#f53dbc")) +
    labs(title="VOC prevalence compared to historic scenario", x="weeks", y=NULL, col=NULL) +
    theme_half_open() +
    theme(legend.position="none")
    # scale_y_continuous(trans="logit")

plot_grid(p1, p2, nrow=1, rel_widths=c(1.2, 1))


ggsave("figs/model-derivation1.pdf", width=10, height=4)
