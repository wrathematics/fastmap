library(rbenchmark)

cols <- c("test", "replications", "elapsed", "relative")
reps <- 5

x = matrix(rnorm(10000*1250), 10000)
benchmark(fastmap:::.fastmap(x), fastmap::fastmap(x), replications=reps, columns=cols)
