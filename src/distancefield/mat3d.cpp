#include "mat3d.h"

#include<stdio.h>

// taken from http://www.dillgroup.ucsf.edu/~bosco/rmsd.c by Bosco K Ho,
// which was, in turn, taken from Numerical Recipes
#define ROTATE(a,i,j,k,l) { g = a[i][j]; \
                            h = a[k][l]; \
                            a[i][j] = g-s*(h+g*tau); \
                            a[k][l] = h+s*(g-h*tau); }

#define SWAP(i,j) { swapEig = eig_val[i];\
                    eig_val[i] = eig_val[j];\
                    eig_val[j] = swapEig;\
                    swapVec = eig_vec[i];\
                    eig_vec[i] = eig_vec[j];\
                    eig_vec[j] = swapVec;}

bool eigen_sym(Mat3d & a, Vec3d & eig_val, Vec3d eig_vec[3], 
               int maxIterations, double epsilon)
{
  int count, k, i, j;
  double tresh, theta, tau, t, sum, s, h, g, c, b[3], z[3];
  int n_rot;

  double swapEig;
  Vec3d swapVec;
                                                                                                                                                             
  double v[3][3];
  double d[3];
                                                                                                                                                             
  /*Initialize v to the identity matrix.*/
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
      v[i][j] = 0.0;
    v[i][i] = 1.0;
  }
                                                                                                                                                             
  /* Initialize   b and d to the diagonal of a */
  for (i=0; i<3; i++)
    b[i] = d[i] = a[i][i];
                                                                                                                                                             
  /* z will       accumulate terms */
  for (i=0; i<3; i++)
    z[i] = 0.0;
                                                                                                                                                             
  n_rot = 0;
                                                                                                                                                             
  /* "maxIterations"   tries */
  for (count=0;   count<maxIterations; count++)
  {
    /* sum off-diagonal     elements */
    sum = 0.0;
    for(i=0; i<2; i++)
    {
      for (j=i+1; j<3; j++)
        sum +=  fabs(a[i][j]);
    }
                                                                                                                                                             
    /* if converged (to machine underflow if epsilon=0.0) */
    if (sum <= epsilon)
    {
      for (int i = 0; i < 3; i++)
      {
        //eig_vec[i] = Vec3d (v[i][0], v[i][1], v[i][2]);
        eig_vec[i] = Vec3d (v[0][i], v[1][i], v[2][i]);
        eig_val[i] = d[i];
      }
      norm(eig_vec[0]);
      norm(eig_vec[1]);
      norm(eig_vec[2]);

      // sort eigenvalues and eigenvectors from largest abs towards smallest
      if (fabs(eig_val[0]) < fabs(eig_val[1]))
        SWAP(0,1);

      if (fabs(eig_val[1]) < fabs(eig_val[2]))
        SWAP(1,2);

      if (fabs(eig_val[0]) < fabs(eig_val[1]))
        SWAP(0,1);

      return true;
    }
                                                                                                                                                             
    /* on 1st three sweeps... */
    if (count < 3)
     tresh = sum * 0.2 / 9.0;
    else
     tresh   = 0.0;
                                                                                                                                                             
    for(i=0; i<2; i++)
    {
      for (j=i+1; j<3; j++)
      {
        g = 100.0 * fabs(a[i][j]);
                                                                                                                                                             
        /*      after four sweeps, skip the     rotation if
         *        the off-diagonal element is small
         */
        if ( count > 3  && fabs(d[i])+g == fabs(d[i]) &&  fabs(d[j])+g == fabs(d[j]) )
        {
          a[i][j] = 0.0;
        }
        else if (fabs(a[i][j]) > tresh)
        {
          h = d[j] - d[i];
                                                                                                                                                             
          if (fabs(h)+g   == fabs(h))
          {
            t =  a[i][j] / h;
          }
          else
          {
            theta = 0.5 * h  / (a[i][j]);
            t = 1.0 / ( fabs(theta) + (double)sqrt(1.0 + theta*theta) );
            if (theta < 0.0)
              t = -t;
          }

          c = 1.0 / (double) sqrt(1 + t*t);
          s = t * c;
          tau = s / (1.0 + c);
          h = t * a[i][j];

          z[i] -= h;
          z[j] += h;
          d[i] -= h;
          d[j] += h;

          a[i][j] = 0.0;

          for (k=0; k<=i-1; k++)
            ROTATE(a, k, i, k, j)
          for (k=i+1; k<=j-1; k++)
            ROTATE(a, i, k, k, j)
          for (k=j+1; k<3; k++)
            ROTATE(a, i, k, j, k)
          for (k=0; k<3; k++)
            ROTATE(v, k, i, k, j)
          ++n_rot;
        }
      }
    }
                                                                                                                                                             
    for(i=0; i<3; i++)
    {
      b[i] += z[i];
      d[i] = b[i];
      z[i] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++)
  {
    //eig_vec[i] = Vec3d (v[i][0], v[i][1], v[i][2]);
    eig_vec[i] = Vec3d (v[0][i], v[1][i], v[2][i]);
    eig_val[i] = d[i];
  }
  norm(eig_vec[0]);
  norm(eig_vec[1]);
  norm(eig_vec[2]);

  // sort eigenvalues and eigenvectors from largest abs towards smallest
  if (fabs(eig_val[0]) < fabs(eig_val[1]))
    SWAP(0,1);

  if (fabs(eig_val[1]) < fabs(eig_val[2]))
    SWAP(1,2);

  if (fabs(eig_val[0]) < fabs(eig_val[1]))
    SWAP(0,1);

  return false;
}
