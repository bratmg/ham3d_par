#include <stdio.h>
#include <math.h>
#include "ham3dtypes.h"
#include "ham3dFunctionDefs.h"
#include <mpi.h>

void computeRHS(GRID *g,SOLN *s,double *l2rho, double *linfrho, int myid)
{
  //
  int i,j,k,m,f,n;
  int f1,f2;
  int is,ie;
  int iface;
  int idv;
  int chainSize;
  double ds[3];
  double leftState[NVAR];
  double rightState[NVAR];
  double leftState0[NVAR];
  double rightState0[NVAR];
  double consVar[NVAR];
  double flux[NVAR];
  double gm1=gamm-1.0;
  double gamma1=gamm;
  double specRadius;
  double faceVel=0.;
  double dsnorm,nynx,nx2ny;
  int node1,node2,node3,node4,leftCell,rightCell,icell;
  double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;  
  double xa,ya,za,xb,yb,zb;
  double pp;
  double th,qt,eps;
  double dscheck[2];
  double ref[3][3];
  int iflag,nbase;
  
  int n1,n2,n3,n4,n5,n6,n7,n8;
  double xc,yc,zc;

  int iface1,iface2,rightCell1,rightCell2,c1,c2,c3;
  int nghost,order,ierr;
  
  order = g->order;


               nghost = 2;
  if(order==5) nghost = 3;

  nbase = (g->nfaces-(g->ncells)/(g->nstrand-1)*g->nstrand)/(g->nstrand-1);
  //
  // zero out residual and spectral radii
  //
  dscheck[0]=dscheck[1]=0;
  for(i=0;i<NVAR*g->ncells;i++) s->r[i]=0.0;
  for(i=0;i<g->ncells;i++) s->sigma[i]=0.0;
  //
  // one loop per chain to evaluate fluxes
  // on all the faces in the chain
  //
  

  for(i=0;i<g->nchains;i++)
    {
      iflag = 0;

      f1=g->faceStartPerChain[i];
      f2=g->faceStartPerChain[i+1];
      m= nghost;

      for(f=f1;f<f2;f++)
	  {
	    iface=g->chainConn[f]; //iface = face index
	    g->cindx[m]=g->faces[8*iface+4]; //cindx = cell index
	    m++;
	  }
      
       
      //
      // add buffer cells to the chain
      //
      if (g->chainConn[f1]==g->chainConn[f2-1])
	  {
     iflag = 0;
	  //
	  // this is a closed chain
	  // make it periodic
	  //
	  f           =  f1+1;
	  iface       =  g->chainConn[f];
	  g->cindx[m] =  g->faces[8*iface+4];
	  m++;
	  chainSize   =  m;
	  m=0;
	  for(f=f2-nghost-1;f<f2-1;f++)
	    {
	      iface       =  g->chainConn[f];
	      g->cindx[m] =  g->faces[8*iface+4];
	      m++;
	    }
	  }
      else
	  {
     iflag = 1;
	  //
	  // this is a open chain
	  // -ve index indicates necessity to create
	  // ghost cells
	  //

     // solid bc 
     if(g->test != 1)
     {
        if(order==5) 
        {
          m--;
          g->cindx[m] = -g->cindx[m];
          m++;
          g->cindx[m] = -g->cindx[m-3];
          m++;
          g->cindx[m] = -g->cindx[m-5];
          chainSize = m+1;
          m = 0;
          g->cindx[m] = -g->cindx[m+5];
          m = 1;
          g->cindx[m] = -g->cindx[m+3];
          m = 2;
          g->cindx[m] = -g->cindx[m+1];
        }
        else
        {
	       m--;
	       g->cindx[m]=-g->cindx[m];
	       m++;
	       g->cindx[m]=-g->cindx[m-3];
	       chainSize=m+1;
	       m=0;
	       g->cindx[m]=-g->cindx[m+3];
	       m=1;
	       g->cindx[m]=-g->cindx[m+1];
        }

        //=========================
        //add for parallization
        //========================
        iface1 = g->chainConn[f1];
        iface2 = g->chainConn[f2-1];
        rightCell1 = g->faces[8*iface1+6];
        rightCell2 = g->faces[8*iface2+6];

        if(rightCell1 == -5)
        {
          if(order==5)
          {
            c3 = g->cindx[5];
            c2 = g->cindx[4];
            c1 = g->cindx[3];

            g->cindx[0]  = -1*(c2);
            g->cindx[1]  = -1*(c2);
            g->cindx[2]  = -1*(c2);
          }
          else
          {
            c2 = g->cindx[3];
            c1 = g->cindx[2];
            g->cindx[0]  = -1*(c2);
            g->cindx[1]  = -1*(c2);
          }
        }

        if(rightCell2 == -5)
        {
          if(order==5)
          {
            c3 = g->cindx[chainSize-6];
            c2 = g->cindx[chainSize-5];
            c1 = g->cindx[chainSize-4];
            g->cindx[chainSize-1]  = -1*(c2);
            g->cindx[chainSize-2]  = -1*(c2);
            g->cindx[chainSize-3]  = -1*(c2);
          }
          else
          {
            c2 = g->cindx[chainSize-4];
            c1 = g->cindx[chainSize-3];
            g->cindx[chainSize-1]  = -1*(c2);
            g->cindx[chainSize-2]  = -1*(c2);
          }
        }

        // prriodic bc at only for strand grid 
        if(g->test ==2)
        {
          iface = g->chainConn[f1];
          
          if(iface>nbase*(g->nstrand-1)-1)
          {
            if(order==5)
            {
              m = chainSize-1;
              g->cindx[2]   = g->cindx[g->nstrand];
              g->cindx[1]   = g->cindx[g->nstrand-1];
              g->cindx[0]   = g->cindx[g->nstrand-2];
              g->cindx[m-2] = g->cindx[3];
              g->cindx[m-1] = g->cindx[4]; 
              g->cindx[m]   = g->cindx[5]; 
            }
            else
            {
              m = chainSize-1;
              g->cindx[1]   = g->cindx[g->nstrand];
              g->cindx[0]   = g->cindx[g->nstrand-1];
              g->cindx[m-1] = g->cindx[2]; 
              g->cindx[m]   = g->cindx[3]; 
            }
          }
        } 
      
      } //not for periodic bc condition         
 
      // periodic bc
      if(g->test==1)
      {
        apply_periodic(&g[0],f1,f2, m);
        chainSize=m+1;
        if(order==5) chainSize=m+2;
      }
    }
    //end:


//=========================================
    //if(myid==0)
    //{
    //if(iflag==1) //open chain 
    //{
//  //  if(rightCell1>=0 || rightCell2>=0)
//  //  {
    //printf("Stopping code\n");
    //printf("rightcell1:%d\n",rightCell1);
    //printf("rightcell2:%d\n",rightCell2);
 
    //printf("chain size=%d\n",chainSize);
    //for(k=0;k<=chainSize-1;k++) 
    //{
    //printf("cell index:%d\n",g->cindx[k]);  
    //}
    ////exit(1);
//  //  }
    //}
    //}
//========================================




//============================================
   for(j=0;j<chainSize;j++)
	{
	  icell=g->cindx[j];
	  //if (icell >=0) 
     if (icell >=0 ||(icell==0&&j==nghost)||(icell==0&&j==chainSize-nghost-1)) // not negative
	  {
	    m=NVAR*icell;
	    for(k=0;k<NVAR;k++)
	    {
		  consVar[k]=s->q[m]; // conservative variables : rho, rho*u, rho*e
		  m++;
	    }
	    g->f[j][0]=consVar[0];
	    g->f[j][1]=consVar[1];
	    g->f[j][2]=consVar[2];
	    g->f[j][3]=consVar[3];
	    g->f[j][4]=consVar[4];
	  }
     if(icell<0||(icell==0&&j==nghost-1)||(icell==0&&j==chainSize-nghost))// icell < 0
	  //else   //icell < 0
	  {
	  //
	  // do ghost cells
	  // based on whether they are on the solid boundary on that
	    if (j < nghost) //0,1
		 {
		  iface=g->chainConn[f1];
		 }
	    else // last
		 {
		  iface=g->chainConn[f2-1];
		 }
	    rightCell=g->faces[8*iface+6];

       // this is for surface bc, so variables are projected
	    if (rightCell==-2)  /* this is a face on solid wall */ 
		 {
		  node1=g->faces[8*iface];
		  node2=g->faces[8*iface+1];
        node3=g->faces[8*iface+2];
        node4=g->faces[8*iface+3];
		  
		  x1=g->x[3*node1];
		  y1=g->x[3*node1+1];
		  z1=g->x[3*node1+2];

		  x2=g->x[3*node2];
		  y2=g->x[3*node2+1];
		  z2=g->x[3*node2+2];

        x3=g->x[3*node3];
		  y3=g->x[3*node3+1];
		  z3=g->x[3*node3+2];

		  x4=g->x[3*node4];
		  y4=g->x[3*node4+1];
		  z4=g->x[3*node4+2];
 
   
        // 3D face normal vector (direction?)
        xa = x1 - x3; xb = x2 - x4;
        ya = y1 - y3; yb = y2 - y4;
        za = z1 - z3; zb = z2 - z4;


        ds[0] = 0.5*(za*yb - ya*zb);
        ds[1] = 0.5*(xa*zb - za*xb);
        ds[2] = 0.5*(ya*xb - xa*yb);


        if(icell>0){printf("Stopping! negative icell at i: %d\n",i);exit(1);}
		  icell=-icell;
		  m=NVAR*icell;
		  for(k=0;k<NVAR;k++)
		  {
		    consVar[k]=s->q[m];
		    m++;
		  }

		  dsnorm=ds[0]*ds[0]+ds[1]*ds[1]+ds[2]*ds[2];

        //make reflection matrix
        ref[0][0] = 1.-2.*ds[0]*ds[0]/dsnorm;
        ref[0][1] = -2.*ds[0]*ds[1]/dsnorm;
        ref[0][2] = -2.*ds[0]*ds[2]/dsnorm;
        ref[1][1] = 1.-2.*ds[1]*ds[1]/dsnorm;
        ref[1][0] = -2.*ds[1]*ds[0]/dsnorm;
        ref[1][2] = -2.*ds[1]*ds[2]/dsnorm;
        ref[2][2] = 1.-2.*ds[2]*ds[2]/dsnorm;
        ref[2][0] = -2.*ds[2]*ds[0]/dsnorm;
        ref[2][1] = -2.*ds[2]*ds[1]/dsnorm;

		  //=======================================================
		  g->f[j][0]  =  consVar[0];
		  g->f[j][1]  =  (consVar[1]*ref[0][0]+consVar[2]*ref[0][1]+consVar[3]*ref[0][2]);
		  g->f[j][2]  =  (consVar[1]*ref[1][0]+consVar[2]*ref[1][1]+consVar[3]*ref[1][2]);
		  g->f[j][3]  =  (consVar[1]*ref[2][0]+consVar[2]*ref[2][1]+consVar[3]*ref[2][2]);
		  g->f[j][4]  =  consVar[4];		  		 
        //========================================================

      }
      //=========================================================
      // add for parallization
      //=========================================================
      else if(rightCell == -5)
      {
           
          icell = -icell; //make positive(second cell index)
          if(j==nghost-2 || j==chainSize-(nghost-1))
          {  
            if(g->idup[icell]==2)
            {
              if(j==nghost-2) c1 = g->cindx[nghost];
              if(j==chainSize-(nghost-1)) c1 = g->cindx[chainSize-(nghost+1)];

              g->f[j][0] = g->dpsil[c1][0];
              g->f[j][1] = g->dpsil[c1][1]; 
              g->f[j][2] = g->dpsil[c1][2]; 
              g->f[j][3] = g->dpsil[c1][3];
              g->f[j][4] = g->dpsil[c1][4];
            }
            else
            { 
              g->f[j][0] = g->psil[icell][0];
              g->f[j][1] = g->psil[icell][1]; 
              g->f[j][2] = g->psil[icell][2]; 
              g->f[j][3] = g->psil[icell][3];
              g->f[j][4] = g->psil[icell][4];
            }

          }
          else if(j==nghost-1 || j==chainSize-(nghost)) //fist cell 
          {
            if(g->idup[icell]==2)
            {
              if(j==nghost-1) c1 = g->cindx[nghost];
              if(j==chainSize-nghost) c1 = g->cindx[chainSize-(nghost+1)];

              g->f[j][0] = g->dpsila[c1][0];
              g->f[j][1] = g->dpsila[c1][1]; 
              g->f[j][2] = g->dpsila[c1][2]; 
              g->f[j][3] = g->dpsila[c1][3];
              g->f[j][4] = g->dpsila[c1][4];
            }
            else
            {
              g->f[j][0] = g->psila[icell][0];
              g->f[j][1] = g->psila[icell][1]; 
              g->f[j][2] = g->psila[icell][2]; 
              g->f[j][3] = g->psila[icell][3];
              g->f[j][4] = g->psila[icell][4];
            }
          }

          if(order==5 && (j==0 || j==chainSize-1))
          {
            //printf("j:%d\n",j);
            if(g->idup[icell]==2)
            {
              if(j==0) c1 = g->cindx[3];
              if(j==chainSize-1) c1 = g->cindx[chainSize-4];

              g->f[j][0] = g->dpsilb[c1][0];
              g->f[j][1] = g->dpsilb[c1][1]; 
              g->f[j][2] = g->dpsilb[c1][2]; 
              g->f[j][3] = g->dpsilb[c1][3];
              g->f[j][4] = g->dpsilb[c1][4];
            }
            else
            { 
              g->f[j][0] = g->psilb[icell][0];
              g->f[j][1] = g->psilb[icell][1]; 
              g->f[j][2] = g->psilb[icell][2]; 
              g->f[j][3] = g->psilb[icell][3];
              g->f[j][4] = g->psilb[icell][4];
            }
          }
		}
	   else // this is for far field b.c, so constant variabls 
		{
        if(g->test==1) 
        {
          printf("Periodic bc has a problem!\n");
          exit(1);
        }
		  g->f[j][0]  =  rinf;
		  g->f[j][1]  =  rinf*s->uinf;
		  g->f[j][2]  =  rinf*s->vinf;
        g->f[j][3]  =  rinf*s->winf;
		  g->f[j][4]  =  pinf/gm1+0.5*rinf*(s->uinf*s->uinf+s->vinf*s->vinf+s->winf*s->winf);
		}
	 } //icell
	 if(g->f[j][0]==0) printf("ERROR::::\n");
	}//chainsize


///////////////////////////////////////////////////////////////////
///////////////end of storing flow varibles ///////////////////////
///////////////////////////////////////////////////////////////////     

   is=nghost-1;
   ie=chainSize-1;
   th=1./3;
   qt=0.25;
   if (g->order==1) qt=0.0;
   eps=1e-10;

   //reconstruction for each chain
   if(order==1 || order==3) 
   {
     muscld(g->f,g->ql,g->qr,g->f2,is,ie,th,qt,eps,chainSize,NVAR);
   }
   if(order==5) 
   {
     weno(g->f,g->ql,g->qr,is,ie,eps,chainSize,NVAR); 
     
     //crweno(g->f,g->ql,g->qr,is,ie,eps,chainSize,NVAR);    
   }
   n=is;
   idv=(g->chainConn[f1]==g->chainConn[f2-1]);

   //
   for(f=f1;f<f2-idv;f++)
	{
	  iface=g->chainConn[f];
	  node1=g->faces[8*iface];
	  node2=g->faces[8*iface+1];
	  node3=g->faces[8*iface+2];
	  node4=g->faces[8*iface+3];

	  leftCell=g->faces[8*iface+4];
	  rightCell=g->faces[8*iface+6];
	  
     //3D case
	  x1=g->x[3*node1];
	  y1=g->x[3*node1+1];
	  z1=g->x[3*node1+2];

	  x2=g->x[3*node2];
	  y2=g->x[3*node2+1];
	  z2=g->x[3*node2+2];

     x3=g->x[3*node3];
	  y3=g->x[3*node3+1];
	  z3=g->x[3*node3+2];

	  x4=g->x[3*node4];
	  y4=g->x[3*node4+1];
	  z4=g->x[3*node4+2];
   
          
     // 3D face normal vector (direction?)
     xa = x3 - x1; xb = x2 - x4;
     ya = y3 - y1; yb = y2 - y4;
     za = z3 - z1; zb = z2 - z4;
     ds[0] = 0.5*(za*yb - ya*zb);
     ds[1] = 0.5*(xa*zb - za*xb);
     ds[2] = 0.5*(ya*xb - xa*yb);

	  for(m=0;m<NVAR;m++)
	  {
	    if (f==f2-idv-1 && idv==0) 
	    {
		  leftState[m]=g->ql[n][m];
		  rightState[m]=g->qr[n+1][m];
		 }
	    else
		 {
		  leftState[m]=g->qr[n+1][m];
		  rightState[m]=g->ql[n][m];
		 }	      
	    leftState0[m]=s->q[NVAR*leftCell+m]; //g->ql[j][m];	            
	   
	    if (rightCell > -1) 
		 {
		  rightState0[m]=s->q[NVAR*rightCell+m]; //g->qr[j+1][m];
		 }  
	  }
	  
     if (rightCell==-1 && g->test!=1) 
     {
	       rightState0[0]=rinf;
	       rightState0[1]=rinf*s->uinf;
          rightState0[2]=rinf*s->vinf;
          rightState0[3]=rinf*s->winf;
          rightState0[4]=pinf/gm1+0.5*rinf*(s->uinf*s->uinf+
          s->vinf*s->vinf+s->winf*s->winf);
     }
	  else if (rightCell==-2 && g->test!=1) 
	  {
       dsnorm=ds[0]*ds[0]+ds[1]*ds[1]+ds[2]*ds[2];
       //make reflection matrix
       ref[0][0] = 1.-2.*ds[0]*ds[0]/dsnorm;
       ref[0][1] = -2.*ds[0]*ds[1]/dsnorm;
       ref[0][2] = -2.*ds[0]*ds[2]/dsnorm;
       ref[1][1] = 1.-2.*ds[1]*ds[1]/dsnorm;
       ref[1][0] = -2.*ds[1]*ds[0]/dsnorm;
       ref[1][2] = -2.*ds[1]*ds[2]/dsnorm;
       ref[2][2] = 1.-2.*ds[2]*ds[2]/dsnorm;
       ref[2][0] = -2.*ds[2]*ds[0]/dsnorm;
       ref[2][1] = -2.*ds[2]*ds[1]/dsnorm;

		 // which fomular is in ?....is it variables? flux?
		 rightState0[0]  =  leftState0[0];
		 rightState0[1]  =  (leftState0[1]*ref[0][0]+leftState0[2]*ref[0][1]+leftState0[3]*ref[0][2]);
		 rightState0[2]  =  (leftState0[1]*ref[1][0]+leftState0[2]*ref[1][1]+leftState0[3]*ref[1][2]);
		 rightState0[3]  =  (leftState0[1]*ref[2][0]+leftState0[2]*ref[2][1]+leftState0[3]*ref[2][2]);
		 rightState0[4]  =  leftState0[4];		  		  

	  }
	  
	  if (rightState[1]!=rightState0[1] && 0)
	  {
	    trace(leftCell);
	    trace(rightCell);
            for(k=0;k<chainSize;k++)
              printf("%d %d %.16e\n",k,g->cindx[k],g->f[k][1]);
	    trace(n);
	    tracef(leftState[1]);
	    tracef(leftState0[1]);
	    tracef(rightState[1]);
	    tracef(rightState0[1]);
	    tracef(s->q[leftCell*NVAR+1])
	    trace(i);
            trace(n);
	    trace(f);
	    exit(0);
	  }
	  //
	  //calculat flux at the face (3D case)

     //fortran 
     flux_roe3d_(ds,leftState,rightState,flux,&specRadius,&gamma1);

	  //
	  m=NVAR*leftCell;
	  for(j=0;j<NVAR;j++)
	  {
	    s->r[m]-=flux[j]; //residual
	    m++;
	  }
	  s->sigma[leftCell]+=specRadius;

	  if (rightCell > -1) 
	  {
	    //printf("rightcell:%d\n",rightCell);
	    m=NVAR*rightCell;
	    for(j=0;j<NVAR;j++)
       {
	     s->r[m] += flux[j]; //residual
        m++;
       }
	    s->sigma[rightCell]+=specRadius;
	  }
	  n++;   
      
	  }  // f1-f2

    } // nstrand
    
  *linfrho=0.;  
  *l2rho=0.;
  for(i=0;i<g->ncells;i++)
    {
      if((*linfrho) < fabs(s->r[5*i]))
      {
       icell=i;
       *linfrho=fabs(s->r[5*i]);
      }
      *l2rho = *l2rho + (s->r[5*i])*(s->r[5*i]);
    }
    *l2rho = sqrt(*l2rho/(g->ncells));
  //tracef(s->sigma[0]);
  //trace(icell);
}
      
      
  


