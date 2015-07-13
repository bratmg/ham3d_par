#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#define nearBodyRadius 3.0

void preprocess(GRID *g, int myid)
{
  int i,m,i1,i2,i3,i4,kk;
  int leftCell,rightCell;
  int leftFaceIndx,rightFaceIndx;
  int f,f1,f2;
  int icell,iface;
  double r1,r2,r3,r4,dis;
  double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,dvedge;
  FILE *fp, *fp1, *fp2, *fp3;
  double volmin,volmax;
  double r_mid[2],s_vec[2];
  double xa,xb,ya,yb,za,zb;
  char   fname[80];  
  int j,k,ff1,ff2,nbase,nbloop,tmp;
  int nsidefaces, nbases;

  int neigb[10000];


  //if (!g->visc) smoothGrid(g,200);
  //find_faces(g->conn,&(g->neig),&(g->faces),&(g->nfaces),g->ncells,4);
  /*
  fp=fopen("QuadData/qedges1.dat","w");
  fprintf(fp,"%d\n",(g->nfaces));
  for(i=0;i<g->nfaces;i++)
    fprintf(fp,"%d %d %d %d %d %d\n",(g->faces[6*i]),
	   (g->faces[6*i+1]),
	   (g->faces[6*i+2]),
	   (g->faces[6*i+4]),
	   (g->faces[6*i+3]),
	   (g->faces[6*i+5]));
  fclose(fp);
  */
  /*
  for(i=0;i<g->nfaces;i++)
    {
      leftCell=g->faces[6*i+2];
      leftFaceIndx=g->faces[6*i+3];
  */  
  //trace(g->nfaces);
  //
  // find boundary faces
  // this code-snippet is not general, its 
  // just for external flow problems where the body is constrained
  // within a radius of nearBodyRadius from origin

  // use periodic bc condition  
  if(g->test==1) periodic_bc(&g[0]); 
  // use surface boundary condition
  if(g->test==0) g->nbfaces = g->ncells/(g->nstrand-1);
  // use infinite wing case 
  if(g->test==2) 
  {
     nbases = g->ncells/(g->nstrand-1); // face on one surface layer
     // nsidefaces = the number of faces only side faces
     nsidefaces = g->nfaces - nbases*g->nstrand;

     g->nbfaces=0;
     for (i=0;i<g->nfaces;i++)
     {
       if(g->faces[8*i+6] == -1 && i < nsidefaces)
       {
         i1=g->faces[8*i];
         i2=g->faces[8*i+1];
         i3=g->faces[8*i+2];
         i4=g->faces[8*i+3];
         //
         x1=g->x[3*i1];
         y1=g->x[3*i1+1];

         x2=g->x[3*i2];
         y2=g->x[3*i2+1];

         x3=g->x[3*i3];
         y3=g->x[3*i3+1];

         x4=g->x[3*i4];
         y4=g->x[3*i4+1];

         r_mid[0] = (x1+x2+x3+x4)/4.;
         r_mid[1] = (y1+y2+y3+y4)/4.;

         dis = sqrt(r_mid[0]*r_mid[0]+r_mid[1]*r_mid[1]);
         if(dis<nearBodyRadius) g->nbfaces++;
        } //if
     } //nbface
  } //test=2

  //
  if(g->test!=1) g->bfaces=(int *) malloc(sizeof(int)*g->nbfaces);
  g->vol=(double *) malloc(sizeof(double)*g->ncells);
  g->neig=(int *) malloc(6*sizeof(int)*g->ncells);
  g->c2f=(int *) malloc(6*sizeof(int)*g->ncells);
  g->c2chain=(int *) malloc(3*sizeof(int)*g->ncells);
  //
  for(i=0;i<g->ncells;i++) g->vol[i]=0.;
  //
  m=0;
  for(i=0;i<g->nfaces;i++)
    {
      leftCell      = g->faces[8*i+4];
      rightCell     = g->faces[8*i+6];
      leftFaceIndx  = g->faces[8*i+5];
      rightFaceIndx = g->faces[8*i+7];
        
      g->c2f[6*leftCell+leftFaceIndx]=i;
      if (rightCell > -1) 
      {
	      g->c2f[6*rightCell+rightFaceIndx]=i;
      }

      //
      // collect neighbors
      //
      g->neig[6*leftCell+leftFaceIndx]=rightCell;
      if (rightCell > -1) 
      {
	     g->neig[6*rightCell+rightFaceIndx]=leftCell;
      }

      //
      // indices of face nodes
      //
      i1=g->faces[8*i];
      i2=g->faces[8*i+1];
      i3=g->faces[8*i+2];
      i4=g->faces[8*i+3];

      //
      x1=g->x[3*i1];
      y1=g->x[3*i1+1];
      z1=g->x[3*i1+2];

      x2=g->x[3*i2];
      y2=g->x[3*i2+1];
      z2=g->x[3*i2+2];

      x3=g->x[3*i3];
      y3=g->x[3*i3+1];
      z3=g->x[3*i3+2];

      x4=g->x[3*i4];
      y4=g->x[3*i4+1];
      z4=g->x[3*i4+2];


      r_mid[0] = (x1+x2+x3+x4)/4.;
      r_mid[1] = (y1+y2+y3+y4)/4.;
      r_mid[2] = (z1+z2+z3+z4)/4.;

      xa = x3-x1;
      xb = x2-x4;
      ya = y3-y1;
      yb = y2-y4;
      za = z3-z1;
      zb = z2-z4;

      s_vec[0]=0.5*(za*yb-ya*zb);
      s_vec[1]=0.5*(xa*zb-za*xb);
      s_vec[2]=0.5*(ya*xb-xa*yb);
     
      dvedge=1./3.*(r_mid[0]*s_vec[0]+r_mid[1]*s_vec[1]+r_mid[2]*s_vec[2]);

      //
      // compute cell volumes
      //
      // hexahedron volume by formula of Gauss
      g->vol[leftCell]+=dvedge;
      
      if (rightCell > -1) g->vol[rightCell]-=dvedge;

      // test=2 : 3D wing
      if(g->test==2)
      { // if a boundary face
        if(g->faces[8*i+6] == -1 && i < nsidefaces)
        {
         dis = sqrt(r_mid[0]*r_mid[0]+r_mid[1]*r_mid[1]);
         if (dis<nearBodyRadius)
          {
            g->bfaces[m]     = i;
            m++;
            g->faces[8*i+6] = -(g->visc+2);
            icell           = g->faces[8*i+4];       // left cell indx
            iface           = g->faces[8*i+5];       // left face 
            g->neig[6*icell+iface] = -(g->visc+2); //??
          }
        }
      }
    }//nface

    //Initial value for volmin and volmax
    volmin=1E15;
    volmax=0.;


   //// data for domain decomp ===================================================

   //fp3 = fopen("./output/conn_total.dat","w");
   //fp2 = fopen("./output/cell_conn.dat","w");
   //fp = fopen("./output/adjncy.dat","w");
   //fp1 = fopen("./output/xadj.dat","w");

   //==========================================================================
   // total cells
   //===========================================================================

   //for(i=0;i<g->ncells;i++)
   //{
   //  fprintf(fp1,"%d\n",kk); //xadj.dat

   //  fprintf(fp3,"%d %d %d %d %d %d\n",g->neig[i*6],g->neig[i*6+1],g->neig[i*6+2],g->neig[i*6+3],g->neig[i*6+4],g->neig[i*6+5]);

   //  for(m=0;m<6;m++)
   //  {
   //     if(g->neig[i*6+m]>=0) 
   //     {
   //       fprintf(fp,"%d\n",g->neig[i*6+m]); //adjncy.dat
   //       if(g->neig[i*6+m]>i)
   //       {
   //       fprintf(fp2,"%d %d\n",i,g->neig[i*6+m]); //cell_conn
   //       }

   //       kk = kk+1;
   //     }
   //  }
   //  if(i==g->ncells-1) fprintf(fp1,"%d\n",kk); //xadj.dat
 
   //}

   ////===============================================================================
   //// only for the surface mesh(nbcell = nbfaces)
   ////=============================================================================
   //kk=0;
   //for(i=0;i<g->nbfaces;i++) 
   //{ 
   //  fprintf(fp1,"%d\n",kk); //xadj.dat
   //  //------------
   //  // make neigb
   //  //------------
   //  k=0;
   //  for(j=0;j<6;j++)
   //  {       
   //    // only for base level cell
   //    if(g->neig[i*6+j]<g->nbfaces && g->neig[i*6+j]!=-2) 
   //    {
   //      neigb[i*4+k]  = g->neig[i*6+j];
   //      k++;         
   //    }
   //  }
   //  if(k!=4) printf("Error on making boundary cell conn.dat\n");
   //  //--------------
   //  //conn_total.dat
   //  //--------------
   //  fprintf(fp3,"%d %d %d %d\n",neigb[i*4],neigb[i*4+1],neigb[i*4+2],neigb[i*4+3]);
   //  for(m=0;m<4;m++)
   //  {
   //     if(neigb[i*4+m]>=0) 
   //     {
   //       fprintf(fp,"%d\n",neigb[i*4+m]); //adjncy.dat
   //       if(neigb[i*4+m]>i)
   //       {
   //       fprintf(fp2,"%d %d\n",i,neigb[i*4+m]); //cell_conn
   //       }
   //       kk = kk+1;
   //     }
   //  }
   //  if(i==g->nbfaces-1) fprintf(fp1,"%d\n",kk); //xadj.dat
   //}

   //fclose(fp);
   //fclose(fp1);
   //fclose(fp2);
   //fclose(fp3);

   //exit(1);

  //==========================================================================


  //for(i=0;i<g->ncells;i++)
  //{
  //printf("c2f(icell:%d): %d %d %d %d %d %d\n",i+1,g->c2f[6*i],g->c2f[6*i+1],g->c2f[6*i+2],g->c2f[6*i+3],g->c2f[6*i+4],g->c2f[6*i+5]);
  //}
  //exit(1);
  //find maximum and minimum v(lume
  for(i=0;i<g->ncells;i++)
  {
   volmin=(volmin < g->vol[i]) ? volmin : g->vol[i]; 
   volmax=(volmax > g->vol[i]) ? volmax : g->vol[i];
  }
  
  if(myid==0)
  {
  tracef(volmin);
  tracef(volmax);
  tracef(g->vol[0]);
  }

  //=============================================================
  // find cell to chain connectivity
  //============================================================
  for(i=0;i<g->ncells;i++)
    { 
      g->c2chain[3*i]=-1;
      g->c2chain[3*i+1]=-1;
      g->c2chain[3*i+2]=-1;
    }
  for(i=0;i<g->nchains;i++)
    {
      f1=g->faceStartPerChain[i];
      f2=g->faceStartPerChain[i+1];

      for(f=f1;f<f2-1;f++)
	   {
	     iface=g->chainConn[f]; //iface = face index that consist chain
	     leftCell=g->faces[8*iface+4]; // left cell

	     if (g->c2chain[3*leftCell] > -1) 
	     {
            if (g->c2chain[3*leftCell+1] >-1) 
            { 
	         g->c2chain[3*leftCell+2]=i;
            }
            else
            {
               g->c2chain[3*leftCell+1]=i;
            }
	     }
	     else
	     {
	         g->c2chain[3*leftCell]=i;
	     }
	  }
   }

   //==================================================================
   // WRITE TO VERIFY (Solid wall face)
   //==================================================================
   if(g->test==2)
   { 
     sprintf(fname,"./output/bface_%d.dat",myid);
     fp = fopen(fname,"w");
     for (i=0;i<g->nbfaces;i++)
     {
       iface = g->bfaces[i];
       i1=g->faces[8*iface];
       i2=g->faces[8*iface+1];
       i3=g->faces[8*iface+2];
       i4=g->faces[8*iface+3];
       fprintf(fp," %f %f %f %f\n",g->x[3*i1],g->x[3*i1+1],g->x[3*i1+2]);
       fprintf(fp," %f %f %f %f\n",g->x[3*i2],g->x[3*i2+1],g->x[3*i2+2]);
       fprintf(fp," %f %f %f %f\n",g->x[3*i3],g->x[3*i3+1],g->x[3*i3+2]);
       fprintf(fp," %f %f %f %f\n",g->x[3*i4],g->x[3*i4+1],g->x[3*i4+2]);

     }
     fclose(fp);
   }

}
		

