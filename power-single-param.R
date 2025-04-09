library(tidyverse)
library(cowplot)

source("src/growth-model.R")
source("src/nopt-single-param.R")

# i.e. "right sided" hypothesis r > 0
# TODO: notice level is reversed from orig. too
nopt_power_right <- function(tmax, q0, r, power=0.8, level=0.05) {
    s2 <- var_single_param(tmax, q0, r)
    Z1 <- qnorm(1-level) # notice one sided this time
    Z2 <- qnorm(power)
    s2 / r^2 * (Z1 + Z2)^2
}

q0 <- 1e-4
r <- 0.11
tmax <- 60

nopt_power_right(tmax, q0, r, power=0.8)
nopt_single_param(tmax, q0, r, error=0.05)

# two population version
# tests against null that both growth rates equal, thus, power gives the probability that
# we correctly conclude the growth rates are different, or the chance of correctly finding 
# evidence against an intrinsic growth advantage

# TODO: think more about how the interpretation to growth advantages would work...dig up some refs
# on this? I know you've seen it mentioned...