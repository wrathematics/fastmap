#' fastmap.nthreads
#' 
#' Returns the number of cores + hyperthreads on the system.  The function
#' respects the environment variable \code{OMP_NUM_THREADS}.
#' 
#' @return
#' An integer; the number of threads.
#' 
#' @export
fastmap.nthreads <- function()
{
  .Call(R_fastmap_nthreads)
}
