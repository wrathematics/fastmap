library(fastmap)

k = 1
x = matrix(1:30, 10)

set.seed(1234)
test = cma(x, k)

set.seed(1234)
truth = fastmap:::.cma(x, k)



stopifnot(all.equal(test, truth))
