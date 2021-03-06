#include <stdio.h>
#include <math.h>
#include "ham3dtypes.h"
#include "ham3dFunctionDefs.h"
void apply_periodic_LHS(GRID *g,int f1, int f2, int m)
{  
  int i,j,k,f,n;
  int is,ie;
  int iface;
  int idv;
  double leftState[NVAR];
  double rightState[NVAR];
  double leftState0[NVAR];
  double rightState0[NVAR];
  int node1,node2,node3,node4,leftCell,rightCell,icell;
  double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;  
  double xa,ya,za,xb,yb,zb;
  double pp;
  double th,qt,eps;
  int iflag,nbase;
  int order, nghost;
  
  order = g->order;
               nghost = 2;
  if(order==5) nghost = 3;

  nbase = (g->nfaces-(g->ncells)/(g->nstrand-1)*g->nstrand)/(g->nstrand-1);
  iface = g->chainConn[f1];

 //strand grid
  if(iface>nbase*(g->nstrand-1)-1)
  {
    g->cindx[1]   = g->cindx[g->nstrand];
    g->cindx[0]   = g->cindx[g->nstrand-1];
    g->cindx[m-1] = g->cindx[2]; 
    g->cindx[m]   = g->cindx[3]; 
    g->ctype[1]   = 1;
    g->ctype[0]   = 1;
    g->ctype[m-1] = 1;
    g->ctype[m]   = 1;
    return;
  } 

   iflag = 0;
   for(j=0;j<g->m1*(g->nstrand-1);j++)
   {
     if(iface == g->bf1[j]) 
     {g->cindx[nghost-1] = g->faces[8*g->bf3[j]+4];
      g->cindx[nghost-2] = g->faces[8*g->bf3_neigf[j][1]+4];
      g->ctype[nghost-1]=1;
      g->ctype[nghost-2]=1;      
      if(order==5) 
      {
       g->cindx[nghost-3] = g->faces[8*g->bf3_neigf[j][2]+4]; 
       g->ctype[nghost-3]=1;
      }
      iflag = 1;
      break;
     }
     if(iface == g->bf3[j]) 
     {g->cindx[nghost-1] = g->faces[8*g->bf1[j]+4];
      g->cindx[nghost-2] = g->faces[8*g->bf1_neigf[j][1]+4];            
      g->ctype[nghost-1]=1;
      g->ctype[nghost-2]=1;

      if(order==5) 
      {
       g->cindx[nghost-3] = g->faces[8*g->bf1_neigf[j][2]+4]; 
       g->ctype[nghost-3]=1;
      }   
      iflag = 1;
      break;
     }
   }

   if(iflag == 0) {
   for(j=0;j<g->m2*(g->nstrand-1);j++)
   {
     if(iface == g->bf2[j]) 
     {g->cindx[nghost-1] = g->faces[8*g->bf4[j]+4];
      g->cindx[nghost-2] = g->faces[8*g->bf4_neigf[j][1]+4];
      g->ctype[nghost-1]=1;
      g->ctype[nghost-2]=1;
      if(order==5) 
      {
       g->cindx[nghost-3] = g->faces[8*g->bf4_neigf[j][2]+4]; 
       g->ctype[nghost-3]=1; 
      }        
      iflag = 1;
      break;
     }
     if(iface == g->bf4[j]) 
     {g->cindx[nghost-1] = g->faces[8*g->bf2[j]+4];
      g->cindx[nghost-2] = g->faces[8*g->bf2_neigf[j][1]+4];
      g->ctype[nghost-1]=1;
      g->ctype[nghost-2]=1;
      if(order==5) 
      {
       g->cindx[nghost-3] = g->faces[8*g->bf2_neigf[j][2]+4]; 
       g->ctype[nghost-3]=1;
      }    
      iflag = 1;
      break;
     }
    }
   }

   //f2-1
   iflag=0;
   iface = g->chainConn[f2-1];
   for(j=0;j<g->m1*(g->nstrand-1);j++)
   {
     if(iface == g->bf1[j]) 
     {g->cindx[m-1] = g->faces[8*g->bf3[j]+4];
      g->cindx[m] = g->faces[8*g->bf3_neigf[j][0]+4];
      g->ctype[m-1]=1;
      g->ctype[m]=1;
      if(order==5) 
      {
       g->cindx[m+1] = g->faces[8*g->bf3_neigf[j][1]+4]; 
       g->ctype[m+1]=1;
      }
      iflag = 1;
      break;
     }
     if(iface == g->bf3[j]) 
     {g->cindx[m-1] = g->faces[8*g->bf1[j]+4];
      g->cindx[m] = g->faces[8*g->bf1_neigf[j][0]+4];
      g->ctype[m-1]=1;
      g->ctype[m]=1;
      if(order==5) 
      {
       g->cindx[m+1] = g->faces[8*g->bf1_neigf[j][1]+4]; 
       g->ctype[m+1]=1;
      }
      iflag = 1;
      break;
     }
    }

    if(iflag == 0){
     for(j=0;j<g->m2*(g->nstrand-1);j++)
     {
     if(iface == g->bf2[j]) 
     {g->cindx[m-1] = g->faces[8*g->bf4[j]+4];
      g->cindx[m] = g->faces[8*g->bf4_neigf[j][0]+4];
      g->ctype[m-1]=1;
      g->ctype[m]=1;
      if(order==5) 
      {
       g->cindx[m+1] = g->faces[8*g->bf4_neigf[j][1]+4]; 
       g->ctype[m+1]=1;
      }
      iflag = 1;
      break;
     }
     if(iface == g->bf4[j]) 
     {g->cindx[m-1] = g->faces[8*g->bf2[j]+4];
      g->cindx[m] = g->faces[8*g->bf2_neigf[j][0]+4];
      g->ctype[m-1]=1;
      g->ctype[m]=1;
      if(order==5) 
      {
       g->cindx[m+1] = g->faces[8*g->bf2_neigf[j][1]+4]; 
       g->ctype[m+1]=1;
      }
      iflag = 1;
      break;
     }
   }
  } 
  if(iflag==0) {
   printf("Periodic bc error!\n");
   exit(1);}

//subroutine end
}
