#include "mrvsmooth.h"
#include "vortex.h"

#define RTOL 1e-8
#define M 10000

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
  //double theta   = a[11];
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
    if(r[ix] > RTOL){
      Z4[ix] = Z3*srat*pow((r[ix]/r3),(-1.0-alpha));
    }else{
      Z4[ix] = 0.0;
    }
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
/*****************************************
* minr (km) 
* maxr (km)
* grid_radius (km)
******************************************/
void mrvs_v_grid(int nx, int ny, double minr, double maxr, double *a, double *grid_radius, double *v){
  double *zeta;//vorticity
  int ix,i;
  double *v2;//velocity on 1D radius
  double *rr;//finer 1d grid for r for integration
  //double M=5000;//resolution of radial integration
  double dr,r;
  int blendorder=9;

  /*--------------------------------*/
  dr = (maxr-minr)/(M-1);//-1 so last rr = maxr
  
  rr = malloc(M*sizeof(double));
  for(ix=0;ix<M;ix++){
    rr[ix]=1000.0*(minr+ix*dr);
  }

  zeta = vor(M,rr,a,blendorder);
  for(ix=0;ix<M;ix++){
    zeta[ix] = rr[ix]*zeta[ix];
  }
  /*--------------------------------*/  
  v2 = malloc(M*sizeof(double));
  trapintegrate(M, zeta, rr, 0.0, &v2[0]);
  for(ix=0;ix<M;ix++){
    if(rr[ix] > RTOL){
      v2[ix] = v2[ix]/rr[ix];
    }else{
      v2[ix] = 0.0;
    }
  }
  /*---------------------------------------------*/
  for(i=0;i<nx*ny;i++){
    r=grid_radius[i]*1000.0;
    v[i]=lininterp(M,rr,v2,r);
  }
  free(v2);
  free(rr);
}
/*****************************************
* minr (km) 
* maxr (km)
* grid_radius (km)
* pc (mb)
******************************************/
void mrvs_p_grid(int nx, int ny, double minr, double maxr, double *a, double pc, double f, double *grid_radius, double *p){
  double rho=1.15;//maybe do the pRT thing here in future
  double *zeta;//vorticity
  int ix,i;
  double *v2;//velocity on 1D radius
  double *p2;//pressure on 1D radius
  double *rr;//finer 1d grid for r for integration
  //double M=5000;//resolution of radial integration
  double dr,r;
  int blendorder=9;

  /*--------------------------------*/
  dr = (maxr-minr)/(M-1);//-1 so last rr = maxr
  
  rr = malloc(M*sizeof(double));
  for(ix=0;ix<M;ix++){
    rr[ix]=1000.0*(minr+ix*dr); //convert km to m
  }

  zeta = vor(M,rr,a,blendorder);
  for(ix=0;ix<M;ix++){
    zeta[ix] = rr[ix]*zeta[ix];
  }
  /*--------------------------------*/  
  v2 = malloc(M*sizeof(double));
  p2 = malloc(M*sizeof(double));
  trapintegrate(M, zeta, rr, 0.0, &v2[0]);
  for(ix=0;ix<M;ix++){
    if(rr[ix] > RTOL){
      v2[ix] = v2[ix]/rr[ix];//finishes integral for velocity from vorticity
      v2[ix] = (v2[ix]*v2[ix]/rr[ix] +f*v2[ix])*rho;//sets up next integral for pressure
    }else{
      v2[ix]=0.0;
    }
  }

  trapintegrate(M, v2, rr, 0.0, &p2[0]);
  /*---------------------------------------------*/
  for(i=0;i<nx*ny;i++){
    r=grid_radius[i]*1000.0;
    p[i]=lininterp(M,rr,p2,r)*0.01 + pc;
  }

  free(v2);
  free(p2);
  free(rr);
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
    if(rr[ix] > RTOL){
      v2[ix] = v2[ix]/rr[ix];
    }else{
      v2[ix] = 0.0;
    }
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
    if(rr[ix] > RTOL){
      v2[ix] = v2[ix]/rr[ix];//finishes integral for velocity from vorticity
      v2[ix] = (v2[ix]*v2[ix]/rr[ix] +f*v2[ix])*rho;//sets up next integral for pressure
    }else{
      v2[ix]=0.0;
    }
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
  int    N=600;// km
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
  //double nm=0.5144;// kt to m/s
  //double vv[3]={34*nm,50*nm,63*nm};
  double alps[11]={0.45, 0.45, 0.45, 0.5,
               0.5,  0.45, 0.45, 0.45,
               0.45, 0.45, 0.45};
  double nt=11;// number of fixes
  int i,j;
  //double pe,rr[4],hrs;
  double rm,vm,alp,Lb,zm,pc;
  double r[N]; // radial grid
  double *vr,*pr;
  double ar[14];
  
  for(i=0;i<N;i++){
    r[i]=1000.0 +i*1000;
  }
  
  for(i=0;i<nt;i++){
    pc = data[i][2]*1e2;
    //pe = data[i][3]*1e2;
    //for(j=5;j<9;j++){
    //  rr[j-5] = data[i][j]*1e3;
    //}
    rm = data[i][8]*1e3;
    vm = data[i][9];
    //hrs = data[i][10];
        
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

int test_grid_out(){
  FILE *fp;
  double slon= 290.0,  elon=310.0;
  double slat=-3.0, elat=15.0;
  double dlat=0.05, dlon=0.05;
  int nx,ny;
  int i,j;
  double *glon,*glat;  //grid lon/lat
  double *ang, *radius;// radians, km
  double *p, *v;//pressure and Holland velocity
  double *vx,*vy;
  double *tau_x, *tau_y;
  double *rv;//Rankine velocity
  double *rp;//Rankine pressure
  double alpha=0.45; //Rankine parameter
  double lat1,lon1;
  double Rm=8.0875793991;//radius of vmax (km)
  double Vm=50.903348087;//vmax (m/s)
  double rho=1.15;//density of air (kg/m^3)
  double pc,pn; //central and outer pressures (mb)
  double cor;//coriolis force;
  double ar[14],Lb=5e3,Zm,rm_m;
  double f0 = 4.39e-5;// Coriolis

  lat1=10.0;//location of current storm
  lon1=300.0;
  
  pc=950.0; //central pressure
  pn=1004.0;//outer pressure

  rm_m=Rm*1000.0;
    Zm = 2.1*Vm/rm_m;
    ar[0]=Zm;
    ar[1]=rm_m;
    ar[2]=Zm;
    ar[3]=rm_m;
    ar[4]=Zm;
    ar[5]=rm_m;
    ar[6]=0.5*(1-alpha);
    ar[7]=alpha;
    ar[8]=Lb;
    ar[9]=Lb;
    ar[10]=Lb;
    ar[11]=300.0;
    ar[12]=0.0;
    ar[13]=0.0;

  nx=(int)((elon-slon)/dlon);
  ny=(int)((elat-slat)/dlat);

  dlon = (elon-slon)/nx; //recalculate now
  dlat = (elat-slat)/ny; //recalculate now

  glon=malloc(nx*sizeof(double));MCK(glon);
  glat=malloc(ny*sizeof(double));MCK(glat);
  
  ang   =malloc(nx*ny*sizeof(double));MCK(ang);
  radius=malloc(nx*ny*sizeof(double));MCK(radius);
  p     =malloc(nx*ny*sizeof(double));MCK(p);
  v     =malloc(nx*ny*sizeof(double));MCK(v);
  vx    =malloc(nx*ny*sizeof(double));MCK(vx);
  vy    =malloc(nx*ny*sizeof(double));MCK(vy);
  tau_x =malloc(nx*ny*sizeof(double));MCK(tau_x);
  tau_y =malloc(nx*ny*sizeof(double));MCK(tau_y);
  rv    =malloc(nx*ny*sizeof(double));MCK(rv);
  rp    =malloc(nx*ny*sizeof(double));MCK(rp);
  
  printf("ang mem =  %fmb\n",nx*ny*sizeof(double)/1024.0/1024.0);

  for(i=0;i<nx;i++){ glon[i] = slon + i*dlon;  }
  for(i=0;i<ny;i++){ glat[i] = slat + i*dlat;  }

  /*-----------------------------------------------------------------------------*/
  angle_radius(ny,nx,lat1,lon1,glat,glon,ang,radius);//angles and radii of every grid point cf to storm location

  for(i=0;i<nx*ny;i++){
    if(ang[i] > 180.0){ printf("angle > 180 %f\n",ang[i]); }
    if(ang[i] < -180.0){ printf("angle < 180  %f\n",ang[i]); }
    if(radius[i] > 100000.0){ printf("radius > 100000m %f\n",radius[i]); }
    if(radius[i] < 0.0){ printf("radius < 0m %f\n",radius[i]); }
  }  
  /*-----------------------------------------------------------------------------*/
  
  holland_p(nx,ny,Vm,Rm,pc,pn,rho,radius,p);
  holland_v(nx,ny,Vm,Rm,pc,pn,rho,radius,v);
  
  mrvs_v_grid(nx,ny, 0.0, 6000.0, ar,         radius,rv);
  mrvs_p_grid(nx,ny, 0.0, 6000.0, ar, pc, f0, radius,rp);

  cor=coriolis(lat1);
  uv(nx,ny, rv, cor,ang, vx,vy);

  stress(nx,ny,rp,vx,vy, tau_x,tau_y);

  //rankine_v(nx,ny,Vm,Rm,alpha,radius,rv);

  fp=fopen("out","w+");
  for(j=0;j<ny;j++){
    for(i=0;i<nx;i++){
      fprintf(fp,"%f %f %f %f %f %f %f %f %f %f\n",glon[i], glat[j], 
              v[i+j*nx], p[i+j*nx], tau_x[i+j*nx], tau_y[i+j*nx], ang[i+j*nx], radius[i+j*nx], rv[i+j*nx], rp[i+j*nx]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

  free(radius);
  free(p);
  free(v);
  free(rp);
  free(rv);
  double dx=0.12;
  ny=6320;
  radius=malloc(ny*sizeof(double));
  p=malloc(ny*sizeof(double));
  v=malloc(ny*sizeof(double));
  rv=malloc(ny*sizeof(double));
  rp=malloc(ny*sizeof(double));
  for(i=0;i<ny;i++){
    radius[i]=i*dx;
  }

  holland_p(1,ny,Vm,Rm,pc,pn,rho,radius,p);
  holland_v(1,ny,Vm,Rm,pc,pn,rho,radius,v);
  
  mrvs_v_grid(1,ny, 0.0, ny*dx, ar,         radius,rv);
  mrvs_p_grid(1,ny, 0.0, ny*dx, ar, pc, f0, radius,rp);

  fp=fopen("out2","w+");
  for(i=0;i<ny;i++){
    fprintf(fp,"%f %f %f %f %f\n",radius[i], v[i], p[i], rv[i], rp[i]);
  }
  
  fclose(fp);

  printf("\n\nTo plot with gnuplot:\n");
  printf("gnuplot\n");
  printf("set pm3d map\n");
  printf("splot \"out\" u 1:2:N\n");
  printf("\n");
  printf("where N is..\n");
  printf("3: Holland radial velocity\n");
  printf("4: pressure\n");
  printf("5: tau_x\n");
  printf("6: tau_y\n");
  printf("7: angle\n");
  printf("8: radial distance from storm\n");
  printf("9: Rankine radial velocity\n");
  /*-----------------------------------------------------------------------------*/
  free(p);
  free(v);
  free(vx);
  free(vy);
  free(tau_x);
  free(tau_y);
  free(rv);
  free(rp);
  free(glon);
  free(glat);
  free(ang);
  free(radius);
  return 0;
}
