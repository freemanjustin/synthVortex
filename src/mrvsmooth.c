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
double getmin(int n, double *r){
  double t=1e15;
  int ix;

  for(ix=0;ix<n;ix++){
    if(r[ix] < t){
      t=r[ix];
    }
  }
  return t;
}
/***********************************************
 *                  / 
 * integrate G(x) = | F(x) dx
 *                  /
 ***********************************************/
void trapintegrate(int N, double *F, double *x, double initial, double *G){


}

/* simple linear interpolation */
double lininterp(int m, double *xm, double *F, double x){
  int i=0;

  if(x<xm[0] || x > xm[m-1]){
    printf("input x out of bounds of input array x=%f (%f - %f)\n",x,xm[0],xm[m-1]);
    exit(EXIT_FAILURE);
  }
  while(x > xm[i] && i < m){
    i++;
  }
  return 0.5*(F[i]-F[i-1]);
}
/***********************************************
 *
 * interpolate from F(xm) to G(xn)
 *
 ***********************************************/
void interpolate(int n, double *xn, int m, double *xm, double *F, double *G){
  int i;

  for(i=0;i<n;i++){
    G[i]=lininterp(m,xm,F,xn[i]);
  }

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

  if(Zamp*Zamp < 1e-10){
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

double* mrvs_v(int n, double *r, double *a, int blendorder){
  double *zeta;//vorticity
  int ix;
  double *v1,*v2;
  double *rr;//finer 1d grid for r for integration
  double m=10;//resolution of integration
  double minr,maxr,dr;

  /*--------------------------------*/
  minr=getmin(n,r);
  maxr=getmax(n,r);
  dr = (maxr-minr)/(m*n-1);//-1 so last rr = maxr
  
  rr = malloc(m*n*sizeof(double));
  for(ix=0;ix<n*m;ix++){
    rr[ix]=minr+ix*dr;
  }

  
  zeta = vor(n*m,rr,a,blendorder);
  for(ix=0;ix<n*m;ix++){
    zeta[ix] = rr[ix]*zeta[ix];
  }
  /*--------------------------------*/
  
  v2 = malloc(m*n*sizeof(double));
  v1 = malloc(  n*sizeof(double));//output to original grid size

  trapintegrate(n*m, zeta, rr, 0.0, &v2[0]);
  for(ix=0;ix<n*m;ix++){
    v2[ix] = v2[ix]/rr[ix];
  }
  interpolate(n,r, n*m,rr, v2, &v1[0]);
  
  free(v2);
  free(zeta);
  free(rr);

  return v1;
}

double* mrvs_p(int n, double *r, double *a, double pc, double f, int blendorder){
  double rho=1.15;//maybe do the pRT thing here in future
  double *zeta;//vorticity
  int ix;
  double *p1,*v2;
  double *rr;//finer 1d grid for r for integration
  double m=10;//resolution of integration
  double minr,maxr,dr;
  
  /*--------------------------------*/
  minr=getmin(n,r);
  maxr=getmax(n,r);
  dr = (maxr-minr)/(m*n-1);//-1 so last rr = maxr
  
  rr = malloc(m*n*sizeof(double));
  for(ix=0;ix<n*m;ix++){
    rr[ix]=minr+ix*dr;
  }

  
  zeta = vor(n*m,rr,a,blendorder);
  for(ix=0;ix<n*m;ix++){
    zeta[ix] = rr[ix]*zeta[ix];
  }
  /*--------------------------------*/

  v2 = malloc(m*n*sizeof(double));
  p1 = malloc(  n*sizeof(double));//output to original grid size

  trapintegrate(n*m, zeta, rr, 0.0, &v2[0]);
  for(ix=0;ix<n*m;ix++){
    v2[ix] = v2[ix]/rr[ix];//finishes integral for velocity from vorticity
    v2[ix] = (v2[ix]*v2[ix]/rr[ix] +f*v2[ix])*rho;//sets up next integral for pressure
  }

  trapintegrate(n*m, v2, rr, 0.0, &p1[0]);
  interpolate(n,r, n*m,rr, v2, &p1[0]);
  for(ix=0;ix<n*m;ix++){
    p1[ix] += pc;
  }

  free(v2);
  free(zeta);
  free(rr);

  return p1;
}
