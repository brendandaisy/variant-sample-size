library(tidyverse)

simple_model <- function(q0, A, num_iter) {
    res <- accumulate(1:num_iter, \(q, t) c((A %*% q) / sum(A %*% q)), .init=q0)
    imap_dfr(res, \(q, i) tibble(t=i-1, strain=names(q0), q=q))
}

q0 <- c(0.99, 0.01)
names(q0) <- c("A", "B")
A <- matrix(c(
    0.55, 0.05,
    0.3, 0.1
), nrow=2)

res <- simple_model(q0, A, 20)

res |> 
    ggplot(aes(t, q, col=strain)) +
    geom_line() +
    geom_point()

Abal <- matrix(c(
    0.475, 0.05,
    0, 0.475
), nrow=2, byrow=TRUE)

res <- map_dfr(seq(0.05, 0.45, length.out=6), \(d) {
    A <- Abal
    A[1, 1] <- A[1, 1] - d
    A[2, 1] <- d
    
    r1 <- simple_model(q0, A, 25) |> 
        mutate(d=d, type="off diagonal")
    
    A[2, 1] <- 0
    A[2, 2] <- A[2, 2] + d
    
    r2 <- simple_model(q0, A, 25) |> 
        mutate(d=d, type="diagonal")
    
    bind_rows(r1, r2)
})

res |> 
    filter(strain == "B") |> 
    ggplot(aes(t, q, col=as.factor(d))) +
    geom_line(linewidth=1.05, alpha=0.7) +
    facet_wrap(~type, nrow=1) +
    theme_bw() +
    theme(legend.position="none")
