#include <stdio.h>
#include "ham3dtypes.h"
#include "ham3dFunctionDefs.h"
#include <stdlib.h>
#define NQ 5
#include <mpi.h>

void readGrid(GRID *g, int myid, int nproc)
{
  FILE *fp;
  char coord[40],conn[40],oface[40],qloop[40],iqloop[40],ncolor[40];
  char c, commu[80], domain[80], strand[50];
  int i,l1,l2,l3,j,ii;
  int workSize;
  double fdummy,tmp;
  int ierr, iproc, c1;
  int dummy;
  int k,kk,ndup;
  int *dup;
// =================================================================
// Initialize dimension
// ================================================================
    if(myid<10)
    {
    sprintf(coord, "./QuadData/domain00%i/coord.dat",myid);
    sprintf(conn,  "./QuadData/domain00%i/conn.dat",myid);
    sprintf(oface, "./QuadData/domain00%i/ofaces.dat",myid);
    sprintf(ncolor,"./QuadData/domain00%i/ncolors.dat",myid);
    sprintf(iqloop,"./QuadData/domain00%i/iqloops.dat",myid);
    sprintf(qloop, "./QuadData/domain00%i/qloops.dat",myid);
    sprintf(commu, "./QuadData/domain00%i/quadCommu00%i.dat",myid,myid);
    sprintf(domain, "./QuadData/domain00%i/subdomain00%i.dat",myid,myid);
    sprintf(strand, "./QuadData/domain00%i/nstrands.dat",myid);
    }
    else if(10<=myid<100)
    {
    sprintf(coord, "./QuadData/domain0%i/coord.dat",myid);
    sprintf(conn,  "./QuadData/domain0%i/conn.dat",myid);
    sprintf(oface, "./QuadData/domain0%i/ofaces.dat",myid);
    sprintf(ncolor,"./QuadData/domain0%i/ncolors.dat",myid);
    sprintf(iqloop,"./QuadData/domain0%i/iqloops.dat",myid);
    sprintf(qloop, "./QuadData/domain0%i/qloops.dat",myid);
    sprintf(commu, "./QuadData/domain0%i/quadCommu0%i.dat",myid,myid);
    sprintf(domain, "./QuadData/domain0%i/subdomain0%i.dat",myid,myid);
    sprintf(strand, "./QuadData/domain0%i/nstrands.dat",myid);
    }
    else
    {
    sprintf(coord, "./QuadData/domain%i/coord.dat",myid);
    sprintf(conn,  "./QuadData/domain%i/conn.dat",myid);
    sprintf(oface, "./QuadData/domain%i/ofaces.dat",myid);
    sprintf(ncolor,"./QuadData/domain%i/ncolors.dat",myid);
    sprintf(iqloop,"./QuadData/domain%i/iqloops.dat",myid);
    sprintf(qloop, "./QuadData/domain%i/qloops.dat",myid);
    sprintf(commu, "./QuadData/domain%i/quadCommu%i.dat",myid,myid);
    sprintf(domain, "./QuadData/domain%i/subdomain%i.dat",myid,myid);
    sprintf(strand, "./QuadData/domain%i/nstrands.dat",myid);
    }
    
  //=============================================================
  // read coordinates
  //=============================================================
  fp=fopen(coord,"r");
  fscanf(fp,"%d",&(g->nnodes));

  g->x=(double *) malloc(sizeof(double)*3*g->nnodes);

  for(i=0;i<g->nnodes;i++)
    fscanf(fp,"%lf %lf %lf\n",&(g->x[3*i]),&(g->x[3*i+1]),&(g->x[3*i+2]));
  fclose(fp);

  //=============================================================
  // read connectivity
  //=============================================================
  fp=fopen(conn,"r");
  fscanf(fp,"%d",&(g->ncells));

  g->conn=(int *) malloc(sizeof(int)*8*g->ncells);

  for(i=0;i<g->ncells;i++)    
    fscanf(fp,"%d %d %d %d %d %d %d %d\n",
    &(g->conn[8*i]),&(g->conn[8*i+1]),
    &(g->conn[8*i+2]),&(g->conn[8*i+3]),
    &(g->conn[8*i+4]),&(g->conn[8*i+5]),
    &(g->conn[8*i+6]),&(g->conn[8*i+7]));
  fclose(fp);


  //===========================================================
  // read faces
  //===========================================================
  fp=fopen(oface,"r");
  fscanf(fp,"%d",&(g->nfaces));
  g->faces=(int *) malloc(sizeof(int)*8*g->nfaces);
  for(i=0;i<g->nfaces;i++)
   {
    fscanf(fp,"%d %d %d %d %d %d %d %d\n",&(g->faces[8*i]),//n1
	   &(g->faces[8*i+1]),//n2
	   &(g->faces[8*i+2]),//n3
	   &(g->faces[8*i+3]),//n4
	   &(g->faces[8*i+4]),//c1
	   &(g->faces[8*i+6]),//c2
	   &(g->faces[8*i+5]),//e1
	   &(g->faces[8*i+7]));//e2

      // boundary cell should be -1 at upperlayers
      // this is due to just add the total number of faces at each layer
      // to the previous layer face index
      if(g->faces[8*i+7]==-1){
         g->faces[8*i+6]=-1;
      }

     //
     // swap if left cell = -1  
     //
     if (g->faces[8*i+4]==-1 || g->faces[8*i+4]==-2 || g->faces[8*i+4]==-5) 
	  {
         swap(g->faces[8*i]   , g->faces[8*i+1]);//n2,n4 
         swap(g->faces[8*i+2] , g->faces[8*i+3]);//n2,n4 

	      swap(g->faces[8*i+4] , g->faces[8*i+6]);//cell index
         swap(g->faces[8*i+5] , g->faces[8*i+7]);//element index  
	  }
    }
  fclose(fp);

  //=========================================================
  // read colors
  //=========================================================
  fp=fopen(ncolor,"r");
  fscanf(fp,"%d",&(g->ncolors));
  g->chainsPerColor=(int *) malloc(sizeof(int)*g->ncolors);
  for(i=0;i<g->ncolors;i++)
    fscanf(fp,"%d\n",&(g->chainsPerColor[i]));
  fclose(fp);


  //========================================================
  // read chain information (chain size)
  //========================================================
  fp=fopen(iqloop,"r");
  fscanf(fp,"%d\n",&(g->nchains));
  g->nchains--;
  g->faceStartPerChain=(int *) malloc(sizeof(int)*((g->nchains+1)));
  g->nmaxchain=0;

  ii=0;
  for(i=0;i<g->nchains+1;i++)
    {      
      fscanf(fp,"%d\n",&(g->faceStartPerChain[ii]));
        
      // this routine is for delete the duplicated one
      // at every layers
      if (i > 0) 
	   {
        if(g->faceStartPerChain[ii]==g->faceStartPerChain[ii-1])
        {    
          ii--;
        }
	  g->nmaxchain=max(g->nmaxchain,g->faceStartPerChain[ii]-g->faceStartPerChain[ii-1]);
	  }
    ii++;
   }
  fclose(fp);
  g->nchains =ii-1;

  //
  if(myid==0) trace(g->nmaxchain);
  //

  //========================================================
  // read quad communication data
  //========================================================


   fp=fopen(domain,"r"); //read subdomain.dat
   fscanf(fp,"%d %d\n",&(g->ncommu),&(dummy)); // same number of send and recv
  // printf("ncommu:%d\n",g->ncommu); // same number of send and recv
   
   //skip character line   
   do { c=fgetc(fp);} while (c!='\n');
   
   g->icommu = (int *) malloc(sizeof(int)*8*g->ncommu);
   g->isend  = (int *) malloc(sizeof(int)*2*g->ncommu);
   g->irecv  = (int *) malloc(sizeof(int)*2*g->ncommu);
   g->irecvconn  = (int *) malloc(sizeof(int)*g->ncommu);
   g->irecvconn2  = (int *) malloc(sizeof(int)*6*g->ncommu);
   g->idup  = (int *) malloc(sizeof(int)*g->ncells);
   g->psil   = (double **) malloc(sizeof(double *)*(g->ncells));
   g->psila   = (double **) malloc(sizeof(double *)*(g->ncells));
   g->dpsil   = (double **) malloc(sizeof(double *)*(g->ncells));
   g->dpsila   = (double **) malloc(sizeof(double *)*(g->ncells));
   g->psilb   = (double **) malloc(sizeof(double *)*(g->ncells));
   g->dpsilb   = (double **) malloc(sizeof(double *)*(g->ncells));
   
   g->psil_deri   = (double **) malloc(sizeof(double *)*(g->ncells));
   g->psila_deri   = (double **) malloc(sizeof(double *)*(g->ncells));
   g->dpsil_deri   = (double **) malloc(sizeof(double *)*(g->ncells));
   g->dpsila_deri   = (double **) malloc(sizeof(double *)*(g->ncells));
   
   g->psilb_deri   = (double **) malloc(sizeof(double *)*(g->ncells));
   g->dpsilb_deri   = (double **) malloc(sizeof(double *)*(g->ncells));
   

   for (i=0;i<g->ncells;i++)
   {
     g->psil[i] = (double *) malloc(sizeof(double)*NVAR);
     g->psila[i] = (double *) malloc(sizeof(double)*NVAR);
     g->dpsil[i] = (double *) malloc(sizeof(double)*NVAR);
     g->dpsila[i] = (double *) malloc(sizeof(double)*NVAR);
   
     g->psilb[i] = (double *) malloc(sizeof(double)*NVAR);
     g->dpsilb[i] = (double *) malloc(sizeof(double)*NVAR);

     g->psil_deri[i] = (double *) malloc(sizeof(double)*NVAR);
     g->psila_deri[i] = (double *) malloc(sizeof(double)*NVAR);
     g->dpsil_deri[i] = (double *) malloc(sizeof(double)*NVAR);
     g->dpsila_deri[i] = (double *) malloc(sizeof(double)*NVAR);
   
     g->psilb_deri[i] = (double *) malloc(sizeof(double)*NVAR);
     g->dpsilb_deri[i] = (double *) malloc(sizeof(double)*NVAR);
   
   }
   
   // initialize 
   for (i=0;i<g->ncells;i++)
   for (j=0;j<NVAR;j++)
   {
     g->psil[i][j] = 0.;
     g->psila[i][j] = 0.;
     g->dpsil[i][j] = 0.;
     g->dpsila[i][j] = 0.;
     g->psilb[i][j] = 0.;
     g->dpsilb[i][j] = 0.;
     
     g->psil_deri[i][j] = 0.;
     g->psila_deri[i][j] = 0.;
     g->dpsil_deri[i][j] = 0.;
     g->dpsila_deri[i][j] = 0.;
     g->psilb_deri[i][j] = 0.;
     g->dpsilb_deri[i][j] = 0.;
   }  


   fscanf(fp,"%d\n",&(g->nadjp));


   g->iadjp  = (int *) malloc(sizeof(int)*g->nadjp);
   g->istp   = (int *) malloc(sizeof(int)*g->nadjp);
   g->ilengp = (int *) malloc(sizeof(int)*g->nadjp);
   
      
   for(i=0;i<g->nadjp;i++)
   fscanf(fp,"%d",&(g->iadjp[i]));
   for(i=0;i<g->nadjp;i++)
   fscanf(fp,"%d",&(g->istp[i]));
   for(i=0;i<g->nadjp;i++)
   fscanf(fp,"%d",&(g->ilengp[i]));



   for(i=0;i<g->ncommu;i++)
     fscanf(fp,"%d %d\n",&(g->isend[2*i]),&(g->isend[2*i+1]));
  
   //skip character line   
   do { c=fgetc(fp);} while (c!='\n');

   fscanf(fp,"%d\n",&(g->nadjp));

   for(i=0;i<g->nadjp;i++)
   fscanf(fp,"%d",&(g->iadjp[i]));
   for(i=0;i<g->nadjp;i++)
   fscanf(fp,"%d",&(g->istp[i]));
   for(i=0;i<g->nadjp;i++)
   fscanf(fp,"%d",&(g->ilengp[i]));


   for(i=0;i<g->ncommu;i++)
     fscanf(fp,"%d %d\n",&(g->irecv[2*i]),&(g->irecv[2*i+1]));

   fclose(fp);

   //find duplicate irecv 
   ndup = 0;
   for(i=0;i<g->ncells;i++) g->idup[i] = 0;
   
   for(i=0;i<g->ncommu;i++)
   {
     c1 = g->irecv[i*2];
     if(g->idup[c1]==1) 
     { 
       ndup = ndup+1;
     }  
     g->idup[c1] = g->idup[c1]+1;
   }
   

   fp=fopen(commu,"r"); 
   
   for(i=0;i<g->ncommu;i++)
     fscanf(fp,"%d %d %d %d %d %d %d %d\n",
     &(g->icommu[8*i]),
     &(g->icommu[8*i+1]),
     &(g->icommu[8*i+2]),
     &(g->icommu[8*i+3]),
     &(g->icommu[8*i+4]),
     &(g->icommu[8*i+5]),
     &(g->icommu[8*i+6]),
     &(g->icommu[8*i+7]));
  
     fclose(fp);


   for(i=0;i<g->ncommu;i++)
   {
     c1  = g->irecv[i*2];
     iproc = g->irecv[i*2+1];
     for(j=0;j<g->ncommu;j++)
     {
        if(c1==g->icommu[8*j+1] && iproc == g->icommu[8*j+6])
        {  
           g->irecvconn[i] = g->icommu[8*j+2];
        }
     }
   }


    // make isend global index
    for(i=0;i<g->ncommu;i++)
    {
      c1 = g->isend[2*i];
      iproc = g->isend[2*i+1];
      for(j=0;j<g->ncommu;j++)
      {
        
     //   printf("c1:%d,g->icommu:%d\n",c1,g->icommu[8*j+4]);
        if(c1==g->icommu[8*j+4] && iproc==g->icommu[8*j+6]) 
        {
          g->isend[2*i] = g->icommu[8*j+7];
        }
      }
    }

    dup = (int *) malloc(sizeof(int)*g->ncommu);
    for(i=0;i<g->ncommu;i++) dup[i] = 0;
    for(i=0;i<g->ncommu;i++)
    {
      c1  = g->isend[2*i]; //global icell cell
      iproc = g->isend[2*i+1]; //processor number
      
      for(j=0;j<g->ncommu;j++)
      {
        if(c1 == g->icommu[8*j+7] && iproc == g->icommu[8*j+6] && dup[j]==0)
        { 
           g->irecvconn2[6*i]   = g->icommu[8*j];
           g->irecvconn2[6*i+1] = g->icommu[8*j+1];
           g->irecvconn2[6*i+2] = g->icommu[8*j+2];
           g->irecvconn2[6*i+3] = g->icommu[8*j+3];
           g->irecvconn2[6*i+4] = g->icommu[8*j+4];
           g->irecvconn2[6*i+5] = g->icommu[8*j+5];
           dup[j] = dup[j]+1;
        }
      }
    }
    free(dup);

    if(myid==0) printf("#ham2d: Finished reading communication data.........\n");


//=========================================================
// Finish reading quad communication data
//========================================================

  workSize=g->nmaxchain+5;
  g->ql=(double **)malloc(sizeof(double *)*(workSize));
  g->qr=(double **)malloc(sizeof(double *)*(workSize));
  g->dql=(double **)malloc(sizeof(double *)*(workSize));
  g->dqr=(double **)malloc(sizeof(double *)*(workSize));
  g->f=(double **) malloc(sizeof(double *)*(workSize));
  g->fv=(double **) malloc(sizeof(double *)*(workSize));
  g->df=(double **) malloc(sizeof(double *)*(workSize));
  g->f2=(double **) malloc(sizeof(double *)*(workSize));
  g->cindx=(int *)malloc(sizeof(int)*(workSize));
  g->ctype=(int *)malloc(sizeof(int)*(workSize));
  g->A=(double ***) malloc(sizeof(double **)*(workSize));
  g->B=(double ***) malloc(sizeof(double **)*(workSize));
  g->C=(double ***) malloc(sizeof(double **)*(workSize));  
  g->F=(double **) malloc(sizeof(double *)*(workSize));
  g->Q=(double **) malloc(sizeof(double *)*(workSize));
  //
  for(i=0;i<workSize;i++)
    {
      g->ql[i]=(double *)malloc(sizeof(double)*NVAR);
      g->qr[i]=(double *)malloc(sizeof(double)*NVAR);
      g->dql[i]=(double *)malloc(sizeof(double)*NVAR);
      g->dqr[i]=(double *)malloc(sizeof(double)*NVAR);
      g->f[i]=(double *)malloc(sizeof(double)*NVAR);
      g->fv[i]=(double *)malloc(sizeof(double)*NVAR);
      g->df[i]=(double *)malloc(sizeof(double)*NVAR);
      g->f2[i]=(double *)malloc(sizeof(double)*NVAR);
      g->A[i]=(double **) malloc(sizeof(double *)*NQ);
      g->B[i]=(double **) malloc(sizeof(double *)*NQ);
      g->C[i]=(double **) malloc(sizeof(double *)*NQ);
      g->F[i]=(double *) malloc(sizeof(double)*NQ);
      g->Q[i]=(double *) malloc(sizeof(double)*NQ);
      for(j=0;j<NQ;j++)
	{
	  g->A[i][j]=(double *)malloc(sizeof(double)*NQ);
	  g->B[i][j]=(double *)malloc(sizeof(double)*NQ);
	  g->C[i][j]=(double *)malloc(sizeof(double)*NQ);
	}
    }
  //
  // read chain information (chain connecitivity)
  //
  fp=fopen(qloop,"r");
  fscanf(fp,"%d",&(g->nchainFaces));
  g->chainConn=(int *)malloc(sizeof(int) * g->nchainFaces);
  for(i=0;i<g->nchainFaces;i++)
    fscanf(fp,"%d",&(g->chainConn[i]));
  fclose(fp);

  fp=fopen(strand,"r");
  fscanf(fp,"%d",&(g->nstrand));
  fclose(fp);

  if(myid==0)
  {
  printf("#ham3d: Finished reading files\n");
  trace(g->nnodes);
  trace(g->ncells);
  trace(g->nfaces);
  trace(g->ncolors);
  trace(g->nchains);
  trace(g->nchainFaces);
  trace(g->nstrand);
  }
}

