#include <stdlib.h>// rand
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MCK(d) {if(d==NULL){printf("Malloc failed %s in %s at line %d\n", __FILE__,__func__, __LINE__); exit(1);}}


// prototypes

double coriolis(double lat);
void uv(int nx, int ny, double *vsf, double f, double *ang, double *u, double *v);
void angle_radius(int ny, int nx, 
                  double lat1,  double lon1,
                  double *glat, double *glon,
                  double *ang,  double *radius);
void stress(int nx, int ny, double *p,double *u,double *v, double *tau_x,double *tau_y);
void rankine_v(int nx, int ny, double Vm, double Rm, double alpha, double *radius, double *v);
double hA(double Rm, double B);
double hB(double Vm, double Rm, double pc, double pn, double rho);
void holland_v(int nx, int ny, double Vm, double Rm, double pc, double pn, double rho, double *radius, double *v);
void holland_p(int nx, int ny, double Vm, double Rm, double pc, double pn, double rho, double *radius, double *p);
int test_vortex();


