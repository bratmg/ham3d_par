#include <stdio.h>
// http://en.wikibooks.org/wiki/Algorithm_Implementation/
//                Linear_Algebra/Tridiagonal_matrix_algorithm
void tri(double *a, double *b, double *c, double *d, int n){

  int i;
  double q; //intermediate quantity

  n--; // since we start from x0 (not x1)
  c[0] /= b[0];
  d[0] /= b[0];

  for (i = 1; i < n; i++) {
    //This quantity is used twice, so calculate it once to
    //avoid some flops and array accesses.
    q = b[i] - a[i]*c[i-1];
    c[i] /= q;
    d[i] = (d[i] - a[i]*d[i-1]) / q;
  }

  d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

  for (i = n; i-- > 0;) {
    d[i] -= c[i]*d[i+1];
  }

}
