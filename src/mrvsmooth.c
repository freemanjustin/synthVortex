#include "mrvsmooth.h"

double getmax(int n, double *r){
  double t=-1e15;
  int ix;

  for(ix=0;ix<n;ix++){
    if(r[ix] > t){
      t=r[ix];
    }
  }
  return t;
}
/******************************************************************************/

double* wt(int n, double *r, double ra, double Lblend, int blendorder){
  int ix;
  double *w;
  double rs;
  
  w=malloc(n*sizeof(double));
  for(ix=0;ix<n;ix++){
    rs = (r[ix]-ra)/Lblend; 
    if (blendorder == 5){
      //# 5th-order blending
      w[ix] = 0.5 - 0.0625*rs*(3.0*pow(rs,4.0) - 10.0*pow(rs,2.0) + 15.0);
    }
    else if (blendorder == 9){
      //# 9th-order blending
      w[ix] = 1.0 - pow(1.0 + rs,5.0) * (128.0 - rs*(325.0 - rs*(345.0 - rs*(175.0 - rs*35.0)))) / 256.0;
    }
  }
  
  return w;
}


  double* vor(int n, double *r, double a[14], int blendorder){
  int ix;
  double *vor;
  double *Z4,*Z5;
  double *w1,*w2,*w3,*w4;

  double Z1      = a[0];
  double r1      = a[1];
  double Z2      = a[2];
  double r2      = a[3];
  double Z3      = a[4];
  double r3      = a[5];
  double srat    = a[6];
  double alpha   = a[7];
  double Lblend1 = a[8];
  double Lblend2 = a[9];
  double Lblend3 = a[10];
  double theta   = a[11];
  double Zamp    = a[12];
  double Zlam    = a[13];
  double Lblend4;
  double r4;
  
  if(isnan(Z2)){
    Z2 = Z1;
    r2 = r1;
  }
  if(isnan(Z3)){
    Z3 = Z2;
    r3 = r2;
  }
  
  Z4 = malloc(n*sizeof(double));
  Z5 = malloc(n*sizeof(double));

  for(ix=0;ix<n;ix++){
    Z4[ix] = Z3*srat*pow((r[ix]/r3),(-1.0-alpha));
  }

  if(Zamp == 0.0){
    for(ix=0;ix<n;ix++){
      Z5[ix] = 0.0;
    }
  } else {
    for(ix=0;ix<n;ix++){
      Z5[ix] = Zamp*sin(2.0*M_PI*r[ix]/Zlam);
    }
  }

  Lblend4 = 1e4;
  r4 = getmax(n,r) - Lblend4;
  
  w1 = wt(n,r,r1,Lblend1,blendorder);
  w2 = wt(n,r,r2,Lblend2,blendorder);
  w3 = wt(n,r,r3,Lblend3,blendorder);
  w4 = wt(n,r,r4,Lblend4,blendorder);

  vor = malloc(n*sizeof(double));

  for(ix=0;ix<n;ix++){
    vor[ix] = w1[ix]*Z1 + (w2[ix] - w1[ix])*Z2 + (w3[ix] - w2[ix])*Z3 + (1.0 - w3[ix])*Z4[ix] + w4[ix]*Z5[ix];
  }

  free(Z4);
  free(Z5);
  free(w1);
  free(w2);
  free(w3);
  free(w4);
  
  return vor;
}
