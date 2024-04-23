library(tidyverse)

q_true <- 1/10
y <- rbinom(10000, 56, q_true)
q_mle <- y / 56

df1 <- map_dfr(c(1, 5, 10, 50), \(n) tibble(n, y=rbinom(1000, n, 1/10)))

ggplot(df1, aes(as.factor(n), y)) +
    # geom_point(col="purple", shape=1, alpha=0.7, size=1.03) +
    geom_boxplot(fill=NA, col="purple") +
    labs(x="number of sequences", y="VOC detected") +
    theme_bw()

ggplot(tibble(q_mle), aes(q_mle)) +
    geom_histogram(aes(y=after_stat(density)), bins=15, fill="gray80", col="gray40") +
    geom_vline(xintercept=q_true, col="purple", linetype="dashed", linewidth=1.06) +
    stat_function(fun=~dnorm(.x, q_true, sqrt(q_true * (1 - q_true) / 50)), col="orange", linewidth=1.06) +
    labs(x="estimated prevalence") +
    theme_bw()

sample_error <- function(n) {
    y <- rbinom(20000, n, 1/100)
    q_mle <- y / n
    return(q_mle)
    # return(quantile(q_mle, c(0.025, 0.975)))
    ci <- quantile(q_mle, c(0.05, 0.95), type=1)
    # return(max(abs(ci - 1/100)))
}

set.seed(1234)
res <- map_dfr(ceiling(seq(1, 200, length.out=200)), \(n) tibble(n, err=sample_error(n)))

ggplot(res, aes(n, err)) +
    # stat_function(fun=~sqrt(qnorm(0.975) * (1/10 * 9/10) / .x), col="orange") +
    geom_point() +
    geom_vline(xintercept=35) + 
    labs(x="number of sequences", y="error for 95% confidence") +
    theme_bw()

hist(sample_error(36))
unique(sample_error(36))

tmp <- sample_error(36)
tmp[which(tmp == max(tmp))]

quantile(sample_error(35), c(0.05, 0.95), type=3)
