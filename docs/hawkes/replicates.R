library(epihawkes)

src <- readLines("hawkes_malaria.R")

set.seed(42)

eval(parse(text = src[1:grep("^true <- list", src)+1]))               # scenario + true params
eval(parse(text = src[grep("^simulate_tree <- function", src):(grep("^tree <- simulate_tree", src)-1)]))
eval(parse(text = src[grep("^fit_epihawkes <- function", src):(grep("^fit <- fit_epihawkes", src)-1)]))

res <- do.call(rbind, parallel::mclapply(1:20, mc.cores = 8, FUN = function(s) {
  set.seed(100 + s)
  tr <- simulate_tree(T_max, true); f <- fit_epihawkes(tr$times, n_starts = 3)
  d <- report(f); d$seed <- s; d$n <- length(tr$times); d
}))

write.csv(res, "replicates.csv", row.names = FALSE)

s <- aggregate(estimate ~ parameter, res, function(x) c(mean = mean(x), sd = sd(x), lo = min(x), hi = max(x)))

print(cbind(parameter = s$parameter, true = c(a, eta, m, sigma, t_peak), round(s$estimate, 3)))
