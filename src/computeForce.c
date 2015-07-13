#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"

void computeForce(GRID *g,SOLN *s)
{
  int i,n;
  double x1,x2,y1,y2,z1,z2;
  int iface,node1,node2,node3,node4,icell;
  double rho,rhou,rhov,rhow,e,p;
  double fac;
  double cs,ss;

  s->cx=s->cy=s->cz=s->cl=s->cd=0;

  for (i=0;i<g->nbfaces;i++)
    {
      iface = g->bfaces[i];
      node1 = g->faces[8*iface];
      node2 = g->faces[8*iface+1];
      node3 = g->faces[8*iface+2];
      node4 = g->faces[8*iface+3];       
       
      x1    = g->x[3*node1];
      y1    = g->x[3*node1+1];
      z1    = g->x[3*node1+2];

      x2    = g->x[3*node2];
      y2    = g->x[3*node2+1];
      z2    = g->x[3*node2+2];

      icell = g->faces[8*iface+4]; //left cell
      //
      rho   = s->q[NVAR*icell];
      rhou  = s->q[NVAR*icell+1];
      rhov  = s->q[NVAR*icell+2];
      rhow  = s->q[NVAR*icell+3];

      e     = s->q[NVAR*icell+4];
      //
      p     = (gamm-1)*(e-0.5*(rhou*rhou+rhov*rhov+rhow*rhow)/rho);
      
      // 2D case
//      s->cy-=p*(x2-x1);
//      s->cx+=p*(y2-y1);

    }

  // 3D case did not be implemented yet.
  fac   = 0.5*rinf*s->mach*s->mach;
  s->cy-= fac;
  s->cx+= fac;

  cs    = cos(s->alpha*deg2rad);
  ss    = sin(s->alpha*deg2rad);
  s->cl = s->cy*cs-s->cx*ss;
  s->cd = s->cy*ss+s->cx*cs;      
}


