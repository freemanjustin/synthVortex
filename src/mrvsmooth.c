#include "mrvsmooth.h"

float getmax(int n, float *r){
  float t=-1e15;
  int ix;

  for(ix=0;ix<n;ix++){
    if(r[ix] > t){
      t=r[ix];
    }
  }
  return t;
}
/******************************************************************************/

float* wt(int n, float *r, float ra, float Lblend, int blendorder){
  int ix;
  float *w;
  float rs;
  
  w=malloc(n*sizeof(float));
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


  float* vor(int n, float *r, float a[14], int blendorder){
  int ix;
  float *vor;
  float *Z4,*Z5;
  float *w1,*w2,*w3,*w4;

  float Z1      = a[0];
  float r1      = a[1];
  float Z2      = a[2];
  float r2      = a[3];
  float Z3      = a[4];
  float r3      = a[5];
  float srat    = a[6];
  float alpha   = a[7];
  float Lblend1 = a[8];
  float Lblend2 = a[9];
  float Lblend3 = a[10];
  float theta   = a[11];
  float Zamp    = a[12];
  float Zlam    = a[13];
  float Lblend4;
  float r4;
  
  if(isnan(Z2)){
    Z2 = Z1;
    r2 = r1;
  }
  if(isnan(Z3)){
    Z3 = Z2;
    r3 = r2;
  }
  
  Z4 = malloc(n*sizeof(float));
  Z5 = malloc(n*sizeof(float));

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

  vor = malloc(n*sizeof(float));

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
