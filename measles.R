

q0 <- 0.5
r <- -0.05
tmax <- 30*4
t0 <- 0

growth_model2(tmax, t0, q0, r)
ceiling(nopt_single_param(tmax, t0, q0, r, error=0.15*abs(r)))
# ceiling(find_sample_size(tmax, t0, q0, r, error=0.15*abs(r), known_idx=2) / 7)
