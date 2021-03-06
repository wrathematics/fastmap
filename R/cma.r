#' cma
#' 
#' Implementation of the coordinate mapping algorithm.
#' 
#' @details
#' The method is threaded via OpenMP. This requires \code{nthreads*n} doubles of
#' additional storage.  If \code{n} is very large relative to \code{m}, this
#' could be very costly; in this case, try using a smaller number of threads by
#' setting the environment variable \code{OMP_NUM_THREADS}.
#' 
#' @param x
#' Input data, preferably as a matrix, or optionally as anything coercible to
#' a matrix.
#' @param k
#' The number of coordinates you want mapped.
#' 
#' @return
#' An \code{m}x\code{k} matrix.
#' 
#' @examples
#' \dontrun{
#' library(fastmap)
#' 
#' x = matrix(1:30, 10)
#' cma(x)
#' }
#' 
#' @references
#' Ostrouchov, G., 2009. A Matrix Computation View of FastMap and RobustMap
#' Dimension Reduction Algorithms. SIAM Journal on Matrix Analysis and
#' Applications, 31(3), pp.1351-1360.
#' 
#' @export
cma = function(x, k=1)
{
  if (!is.matrix(x))
    x = as.matrix(x)
  
  if (!is.double(x))
    storage.mode(x) = "double"
  
  .Call(R_cma, x, as.integer(k))
}
