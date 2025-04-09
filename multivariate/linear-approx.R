library(tidyverse)

growth_model <- function(tmax, t0, q0, r) {
    a <- 1/q0 - 1
    c(rep(0, t0), 1 / (1 + a*exp(-r * 0:(tmax-t0))))
}

simulate_data <- function(tmax, n, t0_true, q0_true, r_true, nsim=10000) {
    q_true <- growth_model(tmax, t0_true, q0_true, r_true)
    
    y <- map(1:nsim, \(rep) rbinom(tmax+1, n, q_true))
    return(list(y=y, n=n, t0=t0_true, q0=q0_true, r=r_true, tmax=tmax, qmax=last(q_true)))
}

q <- growth_model(35, 5, 1e-4, 0.12)

ggplot(tibble(t=seq_along(q)-1, q=q), aes(t, q)) +
    geom_line()

sim <- simulate_data(35, 1000, 10, 1e-4, 0.12, nsim=10000)

prop <- (do.call(rbind, sim$y) / sim$n) |> 
    as_tibble()

prop |> 
    pivot_longer(everything(), values_to="prop") |> 
    mutate(t=parse_number(name)) |> 
    ggplot(aes(t, prop, group=t)) +
    geom_point(position=position_jitter(width=0), alpha=0.1, shape=1)

emp_logit <- function(p, ep) {
    log((p + ep) / (1-p+ep))
}

prop_trans <- emp_logit(do.call(rbind, sim$y) / sim$n, 1e-4/2) |> 
    as_tibble() |> 
    pivot_longer(everything(), values_to="prop") |> 
    mutate(t=parse_number(name))

lm <- lm(prop ~ t, prop_trans)
summary(lm)

ggplot(prop_trans, aes(t, prop)) +
    geom_boxplot(aes(group=t)) +
    geom_smooth(method="lm", col="red", linewidth=1.1)

as.double(summarize(y, across(everything(), var))[1,]) # obv not homoskedastic



map_dbl(1:nrow(ymat), )
