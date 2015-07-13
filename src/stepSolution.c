#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#include <string.h>
#include <stdio.h>
#define deps 1e-10
void stepSolution(char *stepType,GRID *g,SOLN *s,double dt,double *l2norm, double *linfnorm, int myid)
{
 double     coef;
 static int istep=0;
 int i,k,l;
 int nsubiter=20;
 double CFL;

 CFL=g->CFL;

 if (strcmp(stepType,"euler")==0) 
 {
   communication(g,s,myid);
   if(g->visc)
   {
     computeRHSv(g,s,l2norm,linfnorm,myid);
   }
   else
   {
     computeRHS(g,s,l2norm,linfnorm,myid);
   }
   coef=1.0;
   updateSoln(s->q,s->q,s->r,s->sigma,dt,CFL,coef,g->ncells);
  }
 else if (strcmp(stepType,"rk3") == 0){
   
   /* RK step 1 */
   communication(g,s,myid);
   computeRHS(g,s,l2norm,linfnorm,myid);
   coef=0.25;
   updateSoln(s->q,s->qt,s->r,s->sigma,dt,CFL,coef,g->ncells);
   coef=8.0/15;
   updateSoln(s->q,s->q,s->r,s->sigma,dt,CFL,coef,g->ncells);
   
   /* RK step 2 */
   communication(g,s,myid);
   computeRHS(g,s,l2norm,linfnorm,myid);
   coef=5.0/12;
   updateSoln(s->qt,s->q,s->r,s->sigma,dt,CFL,coef,g->ncells);
   
   /* RK step 3 */
   communication(g,s,myid);
   computeRHS(g,s,l2norm,linfnorm,myid);
   coef=3.0/4.0;
   updateSoln(s->qt,s->q,s->r,s->sigma,dt,CFL,coef,g->ncells);
 }
 else if (strcmp(stepType,"ADI")==0){
   // communication
   communication(g,s,myid);
   computeRHSk(g,s,l2norm,linfnorm,myid);
   ADI(g,s,s->cflnum,dt,myid);
 }
 else if (strcmp(stepType,"DADI")==0){
   // communication
   communication(g,s,myid);
   
   computeRHSk(g,s,l2norm,linfnorm,myid);
   DADI(g,s,s->cflnum,dt,myid);
 }
 else if (strcmp(stepType,"Gauss-Seidel") == 0) {
   // communication
   communication(g,s,myid);

   computeRHSkv(g,s,l2norm,myid);
   gaussSeidel(g,s,s->cflnum,dt,myid);
  }
 else if (strcmp(stepType,"line-Gauss-Seidel") == 0) {
   // communication
   communication(g,s,myid);
   if (g->visc) 
     {
       computeRHSkv(g,s,l2norm,myid);
     }
   else
     {       
       computeRHSk(g,s,l2norm,linfnorm,myid);
     }
   lineGaussSeidel(g,s,s->cflnum,dt,myid);
  }
 istep++;
}
