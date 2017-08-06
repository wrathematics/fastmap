#ifndef FASTMAP_DNRM3_
#define FASTMAP_DNRM3_


#include <math.h>

// Modification of DNRM2 for more general Minkowski 
// Original code is:
// This version written on 25-October-1982.
//     Modified on 14-October-1993 to inline the call to DLASSQ.
//     Sven Hammarling, Nag Ltd.
// Modified code by Drew Schmidt, 2015
static inline double dnrm3(const int n, double *restrict x, const double p)
{
  double scale = 0., ssq = 1.;
  double nrm = 0.;
  
  if (n < 1) 
    return 0.;
  else if (n == 1) 
    return fabs(x[0]);
  
  
  for (int i=0; i<n; i++)
  {
    if (x[i] != 0)
    {
      double absxi = fabs(x[i]);
      if (scale < absxi)
      {
        ssq = 1. + ssq*pow(scale/absxi, p);
        scale = absxi;
      }
      else
        ssq += pow(absxi/scale, p);
    }
  }
    
  nrm = scale * pow(ssq, 1.0/p);
  
  return nrm;
}


#endif
