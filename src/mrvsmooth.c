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
 *                               / 
 * cumulatively integrate G(x) = | F(x) dx
 *                              /
 ***********************************************/
void trapintegrate(int N, double *F, double *x, double initial, double *G){
  int i;

  G[0]=initial;
  for(i=1;i<N;i++){
    G[i] = G[i-1] + 0.5*(F[i]+F[i-1])*(x[i]-x[i-1]); 
  }

}

/* simple linear interpolation */
double lininterp(int m, double *xm, double *Fm, double x){
  int i=0;
  double dx,df;

  if(x<xm[0] || x > xm[m-1]){
    printf("input x out of bounds of input array x=%f (%f - %f)\n",x,xm[0],xm[m-1]);
    exit(EXIT_FAILURE);
  }
  while(x > xm[i] && i < m){
    i++;
  }
  dx = (x-xm[i-1]);
  df = (Fm[i]-Fm[i-1])/(xm[i]-xm[i-1]);
      
  return Fm[i-1] + df*dx;
}
/***********************************************
 *
 * interpolate from F(xm) to G(xn)
 *
 ***********************************************/
void interpolate(int n, double *xn, int m, double *xm, double *Fm, double *Gn){
  int i;

  for(i=0;i<n;i++){
    Gn[i]=lininterp(m,xm,Fm,xn[i]);
  }

}
/******************************************************************************/

double* wt(int n, double *r, double ra, double Lblend, int blendorder){
  int ix;
  double *w;
  double rs;
  
  w=malloc(n*sizeof(double));
  for(ix=0;ix<n;ix++){
    w[ix]=0.0;
  }
  
  for(ix=0;ix<n;ix++){
    rs = (r[ix]-ra)/Lblend;
    if(rs*rs < 1.0){
      //printf("%f ",rs);
      if (blendorder == 5){
        //# 5th-order blending
        w[ix] = 0.5 - 0.0625*rs*(3.0*pow(rs,4.0) - 10.0*pow(rs,2.0) + 15.0);
      }
      else if (blendorder == 9){
        //# 9th-order blending
        w[ix] = 1.0 - pow(1.0 + rs,5.0) * (128.0 - rs*(325.0 - rs*(345.0 - rs*(175.0 - rs*35.0)))) / 256.0;
      }
    }
    else if(rs < -1.0){
      w[ix] = 1.0;
    }
    /* if(ix<20 || ix > 3980){ */
    /*   printf("%f ",w[ix]); */
    /* } */
    
  }
  //printf("\n");
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
  double *p1,*p2,*v2;
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
  p2 = malloc(m*n*sizeof(double));
  p1 = malloc(  n*sizeof(double));//output to original grid size

  trapintegrate(n*m, zeta, rr, 0.0, &v2[0]);
  for(ix=0;ix<n*m;ix++){
    v2[ix] = v2[ix]/rr[ix];//finishes integral for velocity from vorticity
    v2[ix] = (v2[ix]*v2[ix]/rr[ix] +f*v2[ix])*rho;//sets up next integral for pressure
  }

  trapintegrate(n*m, v2, rr, 0.0, &p2[0]);
  interpolate(n,r, n*m,rr, p2, &p1[0]);
  for(ix=0;ix<n;ix++){
    p1[ix] += pc;
  }

  free(v2);
  free(p2);
  free(zeta);
  free(rr);

  return p1;
}

int testfn(){
  int    N=400;// km
  int blendorder=9;
  double f0 = 4.39e-5;// Coriolis
  double data[11][11] = {{-13.5, 160.5,  971, 1006, 407, 231.5, 102,    56,    41, 36,    0},
                         {-13.7, 158.6,  968, 1006, 444, 250,   125,    65,    41, 38.6,  6},
                         {-14.1, 156.7,  966, 1006, 482, 259.25,148.25, 65,    37, 41.2, 12},
                         {-14.4, 155,	 956, 1006, 444, 291.5, 148.25, 65,    37, 46.3, 18},
                         {-14.9, 153.2,  940, 1006, 444, 310.25,162,    65,    33, 54,   24},
                         {-15.7, 151.7,  938, 1004, 407, 310,   164.25, 65,    28, 54,   30},
                         {-16.5, 149.7,  937, 1003, 370, 384,   171.5,  69.75, 28, 54,   36},
                         {-17.1, 148.1,  931, 1003, 370, 319.25,176,    74.25, 24, 56.6, 42},
                         {-17.5, 146.8,  930, 1002, 333, 305.5, 166.75, 69.25, 22, 56.6, 48},
                         {-17.96,146.03, 929, 1002, 333, 291.5, 152.75, 55.5,  22, 56.6, 50.5},
                         {-18.5, 145,	 958, 1006, 500, 273,   148.25, 46.25, 28, 43.7, 54}};  // yasi stats
  double nm=0.5144;// kt to m/s
  double vv[3]={34*nm,50*nm,63*nm};
  double alps[11]={0.45, 0.45, 0.45, 0.5,
               0.5,  0.45, 0.45, 0.45,
               0.45, 0.45, 0.45};
  double nt=11;// number of fixes
  int i,j;
  double pe,rr[4],rm,vm,hrs,b,alp,Lb,zm,pc;
  double r[N]; // radial grid
  double *vr,*pr;
  double ar[14];
  
  for(i=0;i<N;i++){
    r[i]=1000.0 +i*1000;
  }
  
  for(i=0;i<nt;i++){
    pc = data[i][2]*1e2;
    pe = data[i][3]*1e2;
    for(j=5;j<9;j++){
      rr[j-5] = data[i][j]*1e3;
    }
    rm = data[i][8]*1e3;
    vm = data[i][9];
    hrs = data[i][10];
        
    //b = bs[i];
    alp = alps[i];

    Lb = 05e3;
    zm = 2.1*vm/rm;
    ar[0]=zm;
    ar[1]=rm;
    ar[2]=zm;
    ar[3]=rm;
    ar[4]=zm;
    ar[5]=rm;
    ar[6]=0.5*(1-alp);
    ar[7]=alp;
    ar[8]=Lb;
    ar[9]=Lb;
    ar[10]=Lb;
    ar[11]=300.0;
    ar[12]=0.0;
    ar[13]=0.0;

    vr = mrvs_v(N,r,ar,blendorder);
    pr = mrvs_p(N,r,ar,pc,f0,blendorder);

    if(i==9){
      for(j=0;j<N;j++){
        printf("%d %f %f\n",j,vr[j],pr[j]);
      }
    }
    free(vr);
    free(pr);
  }
  return 0;
}
