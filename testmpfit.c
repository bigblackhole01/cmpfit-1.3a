/* 
 * MINPACK-1 Least Squares Fitting Library
 *
 * Test routines
 * 
 * These test routines provide examples for users to familiarize 
 * themselves with the mpfit library.  They also provide a baseline
 * test data set for users to be sure that the library is functioning
 * properly on their platform.
 *
 * By default, testmpfit is built by the distribution Makefile.
 *
 * To test the function of the mpfit library, 
 *   1. Build testmpfit   ("make testmpfit")
 *   2. Run testmpfit     ("./testmpfit")
 *   3. Compare results of your run with the distributed file testmpfit.log
 *
 * This file contains several test user functions:
 *   1. linfunc() linear fit function, y = f(x) = a - b*x
 *      - Driver is testlinfit()
 *   2. quadfunc() quadratic polynomial function, y = f(x) = a + b*x + c*x^2
 *      - Driver is testquadfit() - all parameters free
 *      - Driver is testquadfix() - linear parameter fixed
 *   3. gaussfunc() gaussian peak
 *      - Driver is testgaussfit() - all parameters free
 *      - Driver is testgaussfix() - constant & centroid fixed
 *           (this routine demonstrates in comments how to impose parameter limits)
 *   4. main() routine calls all five driver functions
 *
 * Copyright (C) 2003,2006,2009,2010, Craig Markwardt
 *
 */

/* Test routines for mpfit library
   $Id: testmpfit.c,v 1.7 2012/01/23 14:52:14 irby Exp $
*/

#include "mpfit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

double gaussrand()
{
    static double V1, V2, S;
    static int phase = 0;
    double X;
    
    if(phase == 0) {
        do {
            double U1 = (double)rand() / 0x7fffffff;
            double U2 = (double)rand() / 0x7fffffff;
            
            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);
        
        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);
    
    phase = 1 - phase;
    
    return X;
}

/* This is the private data structure which contains the data points
   and their uncertainties */
struct vars_struct {
  double *x;
  double *y;
  double *ey;
};

/* Simple routine to print the fit results */
void printresult(double *x, double *xact, mp_result *result) 
{
  int i;

  if ((x == 0) || (result == 0)) return;
  printf("  CHI-SQUARE = %f    (%d DOF)\n", 
	 result->bestnorm, result->nfunc-result->nfree);
  printf("        NPAR = %d\n", result->npar);
  printf("       NFREE = %d\n", result->nfree);
  printf("     NPEGGED = %d\n", result->npegged);
  printf("     NITER = %d\n", result->niter);
  printf("      NFEV = %d\n", result->nfev);
  printf("\n");
  if (xact) {
    for (i=0; i<result->npar; i++) {
      printf("  P[%d] = %f +/- %f     (ACTUAL %f)\n", 
	     i, x[i], result->xerror[i], xact[i]);
    }
  } else {
    for (i=0; i<result->npar; i++) {
      printf("  P[%d] = %f +/- %f\n", 
	     i, x[i], result->xerror[i]);
    }
  }
    
}

/* 
 * linear fit function
 *
 * m - number of data points
 * n - number of parameters (2)
 * p - array of fit parameters 
 * dy - array of residuals to be returned
 * vars - private data (struct vars_struct *)
 *
 * RETURNS: error code (0 = success)
 */
int linfunc(int m, int n, double *p, double *dy, double **dvec, void *vars)
{
  int i;
  struct vars_struct *v = (struct vars_struct *) vars;
  double *x, *y, *ey, f;

  x = v->x;
  y = v->y;
  ey = v->ey;

  for (i=0; i<m; i++) {
    f = p[0] - p[1]*x[i];     /* Linear fit function; note f = a - b*x */
    dy[i] = (y[i] - f)/ey[i];
  }

  return 0;
}

/* Test harness routine, which contains test data, invokes mpfit() */
int testlinfit()
{
  double x[] = {-1.7237128E+00,1.8712276E+00,-9.6608055E-01,
		-2.8394297E-01,1.3416969E+00,1.3757038E+00,
		-1.3703436E+00,4.2581975E-02,-1.4970151E-01,
		8.2065094E-01};
  double y[] = {1.9000429E-01,6.5807428E+00,1.4582725E+00,
		2.7270851E+00,5.5969253E+00,5.6249280E+00,
		0.787615,3.2599759E+00,2.9771762E+00,
		4.5936475E+00};
  double ey[10];
  /*      y = a - b*x    */
  /*              a    b */
  double p[2] = {1.0, 1.0};           /* Parameter initial conditions */
  double pactual[2] = {3.20, 1.78};   /* Actual values used to make data */
  double perror[2];                   /* Returned parameter errors */      
  int i;
  struct vars_struct v;
  int status;
  mp_result result;

  memset(&result,0,sizeof(result));       /* Zero results structure */
  result.xerror = perror;
  for (i=0; i<10; i++) ey[i] = 0.07;   /* Data errors */           

  v.x = x;
  v.y = y;
  v.ey = ey;

  /* Call fitting function for 10 data points and 2 parameters */
  status = mpfit(linfunc, 10, 2, p, 0, 0, (void *) &v, &result);

  printf("*** testlinfit status = %d\n", status);
  printresult(p, pactual, &result);

  return 0;
}

/* 
 * quadratic fit function
 *
 * m - number of data points
 * n - number of parameters (2)
 * p - array of fit parameters 
 * dy - array of residuals to be returned
 * vars - private data (struct vars_struct *)
 *
 * RETURNS: error code (0 = success)
 */
int quadfunc(int m, int n, double *p, double *dy, double **dvec, void *vars)
{
  int i;
  struct vars_struct *v = (struct vars_struct *) vars;
  double *x, *y, *ey;

  x = v->x;
  y = v->y;
  ey = v->ey;

  /* printf ("quadfunc %f %f %f\n", p[0], p[1], p[2]); */

  for (i=0; i<m; i++) {
    dy[i] = (y[i] - p[0] - p[1]*x[i] - p[2]*x[i]*x[i])/ey[i];
  }

  return 0;
}

/* Test harness routine, which contains test quadratic data, invokes
   mpfit() */
int testquadfit()
{
  double x[] = {-1.7237128E+00,1.8712276E+00,-9.6608055E-01,
		-2.8394297E-01,1.3416969E+00,1.3757038E+00,
		-1.3703436E+00,4.2581975E-02,-1.4970151E-01,
		8.2065094E-01};
  double y[] = {2.3095947E+01,2.6449392E+01,1.0204468E+01,
		5.40507,1.5787588E+01,1.6520903E+01,
		1.5971818E+01,4.7668524E+00,4.9337711E+00,
		8.7348375E+00};
  double ey[10];
  double p[] = {1.0, 1.0, 1.0};        /* Initial conditions */             
  double pactual[] = {4.7, 0.0, 6.2};  /* Actual values used to make data */
  double perror[3];		       /* Returned parameter errors */      
  int i;
  struct vars_struct v;
  int status;
  mp_result result;

  memset(&result,0,sizeof(result));          /* Zero results structure */
  result.xerror = perror;	                                      
  for (i=0; i<10; i++) ey[i] = 0.2;       /* Data errors */           

  v.x = x;
  v.y = y;
  v.ey = ey;

  /* Call fitting function for 10 data points and 3 parameters */
  status = mpfit(quadfunc, 10, 3, p, 0, 0, (void *) &v, &result);

  printf("*** testquadfit status = %d\n", status);
  printresult(p, pactual, &result);

  return 0;
}

/* Test harness routine, which contains test quadratic data;

   Example of how to fix a parameter
*/
int testquadfix()
{
  double x[] = {-1.7237128E+00,1.8712276E+00,-9.6608055E-01,
		-2.8394297E-01,1.3416969E+00,1.3757038E+00,
		-1.3703436E+00,4.2581975E-02,-1.4970151E-01,
		8.2065094E-01};
  double y[] = {2.3095947E+01,2.6449392E+01,1.0204468E+01,
		5.40507,1.5787588E+01,1.6520903E+01,
		1.5971818E+01,4.7668524E+00,4.9337711E+00,
		8.7348375E+00};

  double ey[10];
  double p[] = {1.0, 0.0, 1.0};        /* Initial conditions */             
  double pactual[] = {4.7, 0.0, 6.2};  /* Actual values used to make data */
  double perror[3];		       /* Returned parameter errors */      
  mp_par pars[3];                      /* Parameter constraints */          
  int i;
  struct vars_struct v;
  int status;
  mp_result result;

  memset(&result,0,sizeof(result));       /* Zero results structure */
  result.xerror = perror;

  memset(pars, 0, sizeof(pars));       /* Initialize constraint structure */
  pars[1].fixed = 1;                   /* Fix parameter 1 */

  for (i=0; i<10; i++) ey[i] = 0.2;

  v.x = x;
  v.y = y;
  v.ey = ey;

  /* Call fitting function for 10 data points and 3 parameters (1
     parameter fixed) */
  status = mpfit(quadfunc, 10, 3, p, pars, 0, (void *) &v, &result);

  printf("*** testquadfix status = %d\n", status);
  printresult(p, pactual, &result);

  return 0;
}

/* 
 * gaussian fit function
 *
 * m - number of data points
 * n - number of parameters (4)
 * p - array of fit parameters 
 *     p[0] = constant offset
 *     p[1] = peak y value
 *     p[2] = x centroid position
 *     p[3] = gaussian sigma width
 * dy - array of residuals to be returned
 * vars - private data (struct vars_struct *)
 *
 * RETURNS: error code (0 = success)
 */
int gaussfunc(int m, int n, double *p, double *dy, double **dvec, void *vars)
{
  int i;
  struct vars_struct *v = (struct vars_struct *) vars;
  double *x, *y, *ey;
  double xc, sig1,sig2,incre,rise;

  x = v->x;
  y = v->y;
  ey = v->ey;

    sig1 = 11.7*11.7;
    sig2 = 8.8*8.8;
    
    for (i=0; i<m; i++) {
        xc = x[i]-p[2];
        incre=(1-exp(-xc*xc/sig1));rise=incre*incre*incre*sqrt(incre);
        if(i>=32)
            
            dy[i] = (y[i] - p[1]*rise*exp(-xc*xc/sig2) - p[0])/ey[i];
        else
            dy[i]=(y[i]-p[0])/ey[i];
    }

  return 0;
}

/* Test harness routine, which contains test gaussian-peak data */
int testgaussfit(int num)
{
    int i=0;
    double x[128];
    double y[128];
    double t0=3*((double)rand()/0x7fffffff -0.5)+32;
    
    double N0=0*gaussrand()+250;
    double raise=0;
    double y0=0;
    
    double digi_preci=0.052/28;
    
    FILE*f0;
    char line[50];
    char txtname [50];
    __builtin___sprintf_chk (txtname, 0, __builtin_object_size (txtname, 2 > 1 ? 1 : 0), "cs%d.txt", num);
    printf ("[%s]\n",txtname);
    f0=fopen(txtname,"r");
    while(fgets(line, sizeof line, f0)!=((void *)0)) {
       y[i]=atof(line);
       // y[i]=0;
        i++;
    }
    fclose(f0);
    
   // printf("%d\n",num);
    for (i=0;i<128;i++)
    {
        x[i]=i;
        y0=0.05*((double)rand()/0x7fffffff -0.5)+0.1*(sin(i)+sin(2*i)+sin(3*i));
        if(i>=t0)
        {
            raise=1-exp(-(x[i]-t0)*(x[i]-t0)/(11.7*11.7));
        }
        
       y[i]=y[i]*35+N0*raise*raise*raise*sqrt(raise)*exp(-(x[i]-t0)*(x[i]-t0)/(8.8*8.8));
        // y[i]=y0+N0*raise*raise*raise*sqrt(raise)*exp(-(x[i]-t0)*(x[i]-t0)/(8.8*8.8));
        if(y[i]>=0)
        y[i]=y[i]-fmod(y[i],digi_preci);
        else
        y[i]=y[i]+fmod(y[i],digi_preci);
       
    }
    double ey[128];
    double peddle;
    for(i=0;i<t0;i++)
    {peddle=peddle+y[i];}
    peddle=peddle/t0;
    double p[] = {peddle, 245, 30};
    
    double pactual[] = {0.0, N0, t0};
    double perror[3];
    mp_par pars[3];
    struct vars_struct v;
    int status;
    mp_result result;
    
    __builtin___memset_chk (&result, 0,sizeof(result), __builtin_object_size (&result, 0));
    result.xerror = perror;
    
    __builtin___memset_chk (pars, 0,sizeof(pars), __builtin_object_size (pars, 0));
    
     pars[0].fixed = peddle;

    
    for (i=0; i<128; i++) ey[i] = 0.5;
    v.x = x;
    v.y = y;
    v.ey = ey;
    
    
    
    status = mpfit(gaussfunc, 128, 3, p, pars, 0, (void *) &v, &result);
    
    
    p[1]=p[1]*1809/250;
    p[2]=p[2]-t0;
    
    FILE*fp=fopen("delta2.txt","a");
    fprintf(fp,"%f\t%f\n",p[1],p[2]);fclose(fp);
    return 0;

}


/* Test harness routine, which contains test gaussian-peak data 

   Example of fixing two parameter

   Commented example of how to put boundary constraints
*/

/* Main function which drives the whole thing */
int main(int argc, char *argv[])
{
  int i;
  int niter = 1;
  
    double total_time;
   double start = omp_get_wtime();
    
    srand(time(((void *)0)));
    omp_set_dynamic(0);
    omp_set_num_threads(8);
    #pragma omp parallel
    {
        #pragma omp for
     for(int ii=1; ii<1200; ii++)
        testgaussfit(ii);
        
   
       
    
    printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
    
    }
    
    total_time =omp_get_wtime()-start;
    
    printf("\nTime is: %f\n", total_time);
  exit(0);
}
