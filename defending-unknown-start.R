# analysis to justify our decision to consider the start time unknown
# a good supp fig 1
# 1. a figure illustrating the data generation process, motivating: "when did we start?"
# 2. we could analyse proportional data with 3 models: linear transform, simple logistic, logistic w/ unknown start

# supp fig 2: illustration of patching the logistic with unknown start, with q as func of t0

# remember this is just a proof of concept and you will rerun this with 
# "more relevant" values about 1000 times :)

library(tidyverse)

source("src/growth-model.R")

r <- 0.12 * 7 # work with a weekly example
n <- 100

set.seed(123)
sim <- simulate_data(12, n, 3, 1e-3, r, 10000)
ex_dat <- tibble(prop=sim$y[[1]]/n, t=seq_along(prop)-1)

res_full <- voc_growth_ci_mc(sim, nsim=NULL, r_int=c(0, 1))

ml_q <- tibble(
    t=0:12,
    q1=growth_model(12, 3, 1e-3, 0.84),
    q2=growth_model(12, 2, 1e-3, 0.76),
    q3=growth_model(12, 1, 1e-3, 0.64)
) |> 
    pivot_longer(-t)

# TODO: hear me out: logistic growth so does not matter in the regime we care about for this paper
# exponential growth matches so well it takes work to find where logistic growth matters even here.
# if we ARE gonna do a simplified linear analysis, do it with an exponential model. However this 
# probably won't remove the t0 difficulties so who cares

# TODO: growth rate seems to easy. We should have a second, much harder, but equally motivated
# parameter as an example to estimate throughout the paper. (plenty of examples in lit)
    
ggplot(ex_dat, aes(t, prop)) +
    geom_point() +
    geom_line(aes(y=value, col=name), ml_q) +
    labs(x="weeks", y="VOC proportion")

ggplot(res_full$mle, aes(r, t0)) +
    geom_hex() +
    annotate("point", x=r, y=3, col="tomato2", shape=8)
