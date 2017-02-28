#include "vortex.h"

/*-------------------------------------------------------------------------------------
! Calculation of coriolis parameter f= 2*omega*sin(lat)
!-------------------------------------------------------------------------------------*/

double coriolis(double lat){
  double pi=M_PI;
  double omg=7.292e-5;
  double lat_rad=lat*pi/180.0;
  double coriolis;
  return coriolis = 2*omg*sin(lat_rad);
}

/*-------------------------------------------------------------------------------------
! Calculate U and V for forcing file from VSF
! e.g. vsf = Holland radial velocity.
! f is Coriolis force.
!-------------------------------------------------------------------------------------*/
void uv(int nx, int ny, double *vsf, double f, double *ang, double *u, double *v){
  int i,j;
  double f_sign=1.0;
  double pi=M_PI;
  
  if(f<0) f_sign=-1.0;

  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
      u[i+j*nx] = -f_sign*sin(pi/2 - ang[i+j*nx])*vsf[i+j*nx];
      v[i+j*nx] =  f_sign*cos(pi/2 - ang[i+j*nx])*vsf[i+j*nx];
    }
  }

}
/*-------------------------------------------------------------------------------------
! distance field relative to TC centre
! spherical trig great arc seperation and bearing on a sphere.
! output units: radians, kilometres
!-------------------------------------------------------------------------------------*/
void angle_radius(int ny, int nx, 
                  double lat1,  double lon1,
                  double *glat, double *glon,
                  double *ang,  double *radius){
  int i,j;
  double radius_earth=6371.0;
  double pi=M_PI;
  double lat1_rad, lon1_rad;
  double lat2_rad, lon2_rad;
  double dlon_rad, dlat_rad;

  lat1_rad = lat1*pi/180.0;
  lon1_rad = lon1*pi/180.0;

  for(j=0;j<ny;j++){
    lat2_rad = glat[j]*(pi/180.0);
    for(i=0;i<nx;i++){
      lon2_rad = glon[i]*(pi/180.0);
      dlon_rad = lon2_rad-lon1_rad;
      ang[i+j*nx]    = atan2((cos(lat2_rad)*sin(dlon_rad)),
                             (sin(lat2_rad)*cos(lat1_rad)-cos(lat2_rad)*sin(lat1_rad)*cos(dlon_rad)));
      radius[i+j*nx] = radius_earth*acos(sin(lat2_rad)*sin(lat1_rad)+cos(lat2_rad)*cos(lat1_rad)*cos(dlon_rad));
    }
  }

}
/*-------------------------------------------------------------------------------------
! wind stress calculation using bulk formula given by Large and Pond (1981)
! except capped in accordance with Powell et al (2003)
!-------------------------------------------------------------------------------------*/
void stress(int nx, int ny, double *p,double *u,double *v, double *tau_x,double *tau_y){
  double rd = 287.04;
  double t =  273.15 + 30.0;
  int i,j;
  double cd,rho,vv;

  for(j=0;j<ny;j++){
    for(i=0;i<nx;i++){
      rho = (p[i+j*nx]*100)/(rd*t);
      vv = sqrt(u[i+j*nx]*u[i+j*nx] + v[i+j*nx]*v[i+j*nx]);
      
      if(vv < 10.92){
        cd = 1.2e-3;
      }else if(vv < 23.23){
        cd = (0.49 + 0.065*vv)*1e-3;
      }else{
        cd = 2.0e-3;
      }

      tau_x[i+j*nx] = rho*cd*u[i+j*nx]*vv;
      tau_y[i+j*nx] = rho*cd*v[i+j*nx]*vv;
    }
  }
}
/*-------------------------------------------------------------------------------------
! Output:
! Rankine vortex radial velocity (m/s)
! Inputs:
! Vm: maximum radial velocity (m/s)
! Rm: maximum velocity radial distance (km)
! alpha: shape parameter
! radius: radial distances from storm for all grid nodes (km)
!-------------------------------------------------------------------------------------*/
void rankine_v(int nx, int ny, double Vm, double Rm, double alpha, double *radius, double *v){
  int i;
  double r;
  
  for(i=0;i<nx*ny;i++){
    r=radius[i];
    if(r<Rm){
      v[i]=r*Vm/Rm;
    }else{
      v[i]=Vm*pow(r/Rm,-alpha);
    }
  }
}
/*-------------------------------------------------------------------------------------
! Input Rm in kilometres.
!-------------------------------------------------------------------------------------*/
double hA(double Rm, double B){
  return pow(Rm,B);
}
/*-------------------------------------------------------------------------------------
! Inputs:
! Rm in kilometres.
! Vm in kmh
!-------------------------------------------------------------------------------------*/
double hB(double Vm, double Rm, double pc, double pn, double rho){
  return Vm*Vm*rho*exp(1.0)*0.01/(pn-pc);
}
/*-------------------------------------------------------------------------------------
! Output:
! Holland vortex radial velocity (m/s)
! Inputs:
! Vm: maximum radial velocity (m/s)
! Rm: maximum velocity radial distance (km)
! pc: central pressure (mb)
! pn: outer pressure   (mb)
! rho: density (kg/m^3)
! radius: radial distances from storm for all grid nodes (km)
!
! NOTE: could use rho as defined in stress func here
!-------------------------------------------------------------------------------------*/

#define RTOL 1e-6

void holland_v(int nx, int ny, double Vm, double Rm, double pc, double pn, double rho, double *radius, double *v){
    double A,B,V;
    double powr,r;
    int i;
    
    B=hB(Vm,Rm,pc,pn,rho);
    A=hA(Rm,B);
    
    printf("velocity A=%f B=%f\n",A,B);
    
    for(i=0;i<nx*ny;i++){
        r=radius[i];
        if(r>RTOL){
            powr=pow(r,B);
            V=100.0*A*B*(pn-pc)*exp(-A/powr)/(rho*powr);
            V=sqrt(V);
        }else{
            V=0.0;
        }
        v[i]=V;
    }
    
}

/*-------------------------------------------------------------------------------------
! Output:
! Holland vortex pressure (mb)
! Inputs:
! Vm: maximum radial velocity (m/s)
! Rm: maximum velocity radial distance (km)
! pc: central pressure (mb)
! pn: outer pressure   (mb)
! rho: density (kg/m^3)
! radius: radial distances from storm for all grid nodes (km)
!
! NOTE: could use rho as defined in stress func here
!-------------------------------------------------------------------------------------*/
void holland_p(int nx, int ny, double Vm, double Rm, double pc, double pn, double rho, double *radius, double *p){
  double A,B,V;
  double powr;
  int i;

  B=hB(Vm,Rm,pc,pn,rho);
  A=hA(Rm,B);
  printf("pressure A=%f B=%f\n",A,B);

  for(i=0;i<nx*ny;i++){
    powr = pow(radius[i],B);
    p[i] = pc +(pn-pc)*exp(-A/powr);
  }
}

int test_vortex(){
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
  double alpha=0.4; //Rankine parameter
  double lat1,lon1;
  double Rm=31.3;//radius of vmax (km)
  double Vm=56.0;//vmax (m/s)
  double rho=1.15;//density of air (kg/m^3)
  double pc,pn; //central and outer pressures (mb)
  double cor;//coriolis force;

  lat1=10.0;//location of current storm
  lon1=300.0;
  
  pc=930.0; //central pressure
  pn=1004.0;//outer pressure

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

  cor=coriolis(lat1);
  uv(nx,ny,v,cor,ang, vx,vy);

  stress(nx,ny,p,vx,vy, tau_x,tau_y);

  rankine_v(nx,ny,Vm,Rm,alpha,radius,rv);

  fp=fopen("out","w+");
  for(j=0;j<ny;j++){
    for(i=0;i<nx;i++){
      fprintf(fp,"%f %f %f %f %f %f %f %f %f\n",glon[i], glat[j], v[i+j*nx], p[i+j*nx], tau_x[i+j*nx], tau_y[i+j*nx], ang[i+j*nx], radius[i+j*nx], rv[i+j*nx]);
    }
    fprintf(fp,"\n");
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
  free(glon);
  free(glat);
  free(ang);
  free(radius);

  return 0;
}

