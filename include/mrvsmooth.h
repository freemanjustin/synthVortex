#ifndef mrvsmooth_h
#define mrvsmooth_h
#define mrvsmooth_h

#include <stdlib.h>// rand
#include <stdio.h>
#include <string.h>
#include <math.h>

double getmax(int n, double *r);
double getmin(int n, double *r);
double* wt(int n, double *r, double ra, double Lblend, int blendorder);
double* vor(int n, double *r, double a[14], int blendorder);
void mrvs_v_grid(int nx, int ny, double minr, double maxr, double *a, double *grid_radius, double *v);
void mrvs_p_grid(int nx, int ny, double minr, double maxr, double *a, double pc, double f, double *grid_radius, double *p);
double* mrvs_v(int n, double *r, double *a, int blendorder);
double* mrvs_p(int n, double *r, double *a, double pc, double f, int blendorder);
void trapintegrate(int N, double *F, double *x, double initial, double *G);
void interpolate(int n, double *xn, int m, double *xm, double *F, double *G);
int testfn();
int test_grid_out();
#endif
