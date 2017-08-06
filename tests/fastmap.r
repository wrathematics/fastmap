library(fastmap)

x = matrix(1:30, 10)

set.seed(1234)
test = fastmap(x)

set.seed(1234)
truth = fastmap:::.fastmap(x)



stopifnot(all.equal(test$a, truth$a))
stopifnot(all.equal(test$b, truth$b))
