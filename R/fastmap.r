#' fastmap
#' 
#' Implementation of the fastmap algorithm.
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
#' 
#' @return
#' A list of points \code{a} and \code{b}, the latter being the most distant
#' from the former.
#' 
#' @examples
#' \dontrun{
#' library(fastmap)
#' 
#' x = matrix(1:30, 10)
#' fastmap(x)
#' }
#' 
#' @references
#' Ostrouchov, G., 2009. A Matrix Computation View of FastMap and RobustMap
#' Dimension Reduction Algorithms. SIAM Journal on Matrix Analysis and
#' Applications, 31(3), pp.1351-1360.
#' 
#' Faloutsos, C. and Lin, K.I., 1995. FastMap: A fast algorithm for indexing,
#' data-mining and visualization of traditional and multimedia datasets (Vol.
#' 24, No. 2, pp. 163-174). ACM.
#' 
#' @export
fastmap = function(x)
{
  if (!is.matrix(x))
    x = as.matrix(x)
  
  if (!is.double(x))
    storage.mode(x) = "double"
  
  .Call(R_fastmap, x)
}
