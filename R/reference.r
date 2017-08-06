### Reference implementations. Fairly direct translations from the paper

.alldist = function(d, x) sapply(1:nrow(x), function(i) sum((d-x[i, ])^2))

.fastmap = function(x)
{
  d_i = sample(nrow(x), size=1)
  d = x[d_i, ]
  
  dists = .alldist(d, x)
  a_i = which.max(dists)
  a = x[a_i, ]
  
  dists = .alldist(a, x)
  b_i = which.max(dists)
  b = x[b_i, ]
  
  list(a=a, b=b)
}

.cma = function(x, k)
{
  n = ncol(x)
  
  for (i in 1:k)
  {
    x.i = x[, i:n, drop=FALSE]
    ab = .fastmap(x.i)
    
    x.i = sweep(x.i, MARGIN=2, STATS=ab$a, FUN="-")
    
    w = matrix(ab$b - ab$a)
    v = w
    v[1] = w[1] + sign(w[1]) * norm(w, "2")
    v = v / norm(v, "2")
    y = x.i %*% v
    x[, i:n] = x.i - tcrossprod(2*y, v)
  }
  
  x[, 1:k, drop=FALSE]
}
