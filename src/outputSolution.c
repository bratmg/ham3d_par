#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#include <mpi.h>
void outputSolution(GRID *g,SOLN *s,int nn,int myid, int nproc)
{

  int i,j,n;
  int iface,node1,node2,node3,node4;
  int icell;
  double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,rho,rhou,rhov,rhow,e,pp,cp;
  FILE *fp,*fp1;
  char fname[80];

  if        (nn < 10   ) {sprintf(fname,"./output/vol00%d_%d.dat",nn,myid);}
  else if   (nn < 100  ) {sprintf(fname,"./output/vol0%d_%d.dat",nn,myid);}
  else if   (nn < 1000 ) {sprintf(fname,"./output/vol%d_%d.dat",nn,myid);}
  
  
  fp = fopen(fname,"w");

  //default format
  if(s->outform==0 || g->test==2) //3D case only can this NOW
  {
    fprintf(fp,"VARIABLES = \"X\",\"Y\",\"Z\",\"RHO\",\"RHOU\",\"RHOV\",\"RHOW\",\"E\"\n");
    fprintf(fp,"ZONE ZONETYPE=FEBRICK N= %d E= %d DATAPACKING=BLOCK\n",g->nnodes,g->ncells); 
    fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=NODAL, 4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED, 7=CELLCENTERED, 8=CELLCENTERED)\n");
    for(i=0;i<g->nnodes;i++)
      fprintf(fp,"%f\n",g->x[3*i]);
    fprintf(fp,"\n");
    for(i=0;i<g->nnodes;i++)
      fprintf(fp,"%f\n",g->x[3*i+1]);
    fprintf(fp,"\n");
    for(i=0;i<g->nnodes;i++)
      fprintf(fp,"%f\n",g->x[3*i+2]);
    fprintf(fp,"\n");
    for(n=0;n<NVAR;n++)
      for(i=0;i<g->ncells;i++)
        {
         fprintf(fp,"%f\n",s->q[5*i+n]);
        }
    fprintf(fp,"\n");
    for(i=0;i<g->ncells;i++)
      fprintf(fp,"%d %d %d %d %d %d %d %d\n",g->conn[8*i],g->conn[8*i+1],g->conn[8*i+2],g->conn[8*i+3],g->conn[8*i+4],g->conn[8*i+5],g->conn[8*i+6],g->conn[8*i+7]);

    fclose(fp);
  }


  if(s->outform==1 && g->test!=2) //3D cannot this now!
  {
    fprintf(fp,"VARIABLES = \"X\",\"Y\",\"Z\",\"RHO\",\"U\",\"V\",\"W\",\"CP\"\n");
    fprintf(fp,"ZONE ZONETYPE=FEBRICK N= %d E= %d DATAPACKING=BLOCK\n",g->nnodes,g->ncells); 
    fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=NODAL, 4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED, 7=CELLCENTERED, 8=CELLCENTERED)\n");
    for(i=0;i<g->nnodes;i++)
      fprintf(fp,"%f\n",g->x[3*i]);
    fprintf(fp,"\n");
    for(i=0;i<g->nnodes;i++)
      fprintf(fp,"%f\n",g->x[3*i+1]);
    fprintf(fp,"\n");
    for(i=0;i<g->nnodes;i++)
      fprintf(fp,"%f\n",g->x[3*i+2]);
    fprintf(fp,"\n");

    //rho
    for(i=0;i<g->ncells;i++)
    {
      rho = s->q[5*i];
      fprintf(fp,"%f\n",rho);
    }
    // u
    for(i=0;i<g->ncells;i++)
    {
      rho = s->q[5*i];
      rhou = s->q[5*i+1];
      fprintf(fp,"%f\n",rhou/rho);
    }
    // v
    for(i=0;i<g->ncells;i++)
    {
      rho = s->q[5*i];
      rhov = s->q[5*i+2];
      fprintf(fp,"%f\n",rhov/rho);
    }
    // w
    for(i=0;i<g->ncells;i++)
    {
      rho = s->q[5*i];
      rhow = s->q[5*i+3];
      fprintf(fp,"%f\n",rhow/rho);
    }
    // cp
    for(i=0;i<g->ncells;i++)
    {
      rho = s->q[5*i];
      rhou = s->q[5*i+1];
      rhov = s->q[5*i+2];
      rhow = s->q[5*i+3];
      e = s->q[5*i+4];
      pp=(gamm-1)*(e-0.5*(rhou*rhou+rhov*rhov+rhow*rhow)/rho);
      cp=(pp-pinf)/(0.5*s->mach*s->mach); 
      fprintf(fp,"%f\n",cp);
    }

    fprintf(fp,"\n");
    for(i=0;i<g->ncells;i++)
    {
      fprintf(fp,"%d %d %d %d %d %d %d %d\n",g->conn[8*i],g->conn[8*i+1],
        g->conn[8*i+2],g->conn[8*i+3],g->conn[8*i+4],
        g->conn[8*i+5],g->conn[8*i+6],g->conn[8*i+7]);
   }
    fclose(fp);
  }
 
//======================================================================
  // output for surface boundary
  if        (nn < 10   ) {sprintf(fname,"./output/surf00%d_%d.dat",nn,myid);}
  else if   (nn < 100  ) {sprintf(fname,"./output/surf0%d_%d.dat",nn,myid);}
  else if   (nn < 1000 ) {sprintf(fname,"./output/surf%d_%d.dat",nn,myid);}
  fp = fopen(fname,"w");
 
  //default
  if(s->outform==0 && g->test!=2)
  {
  fprintf(fp,"VARIABLES = \"X\",\"Y\",\"Z\",\"RHO\",\"RHOU\",\"RHOV\",\"RHOW\",\"E\"\n");
  fprintf(fp,"ZONE ZONETYPE=FEQUADRILATERAL N= %d E= %d DATAPACKING=BLOCK\n",g->nnodes/g->nstrand,g->nbfaces); 
  fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=NODAL, 4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED, 7=CELLCENTERED, 8=CELLCENTERED)\n");
  for(i=0;i<g->nnodes/g->nstrand;i++)
    fprintf(fp,"%f\n",g->x[3*i]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes/g->nstrand;i++)
    fprintf(fp,"%f\n",g->x[3*i+1]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes/g->nstrand;i++)
    fprintf(fp,"%f\n",g->x[3*i+2]);
  fprintf(fp,"\n");
  for(n=0;n<NVAR;n++)
    for(i=0;i<g->nbfaces;i++)
      {
       fprintf(fp,"%f\n",s->q[5*i+n]);
      }
  fprintf(fp,"\n");
  for(i=0;i<g->nbfaces;i++)
    fprintf(fp,"%d %d %d %d\n",g->conn[8*i],g->conn[8*i+1],g->conn[8*i+2],g->conn[8*i+3]);

  fclose(fp);
  }
  
  //new
  if(s->outform==1 && g->test!=2)
  {
  fprintf(fp,"VARIABLES = \"X\",\"Y\",\"Z\",\"RHO\",\"U\",\"V\",\"W\",\"CP\"\n");
  fprintf(fp,"ZONE ZONETYPE=FEQUADRILATERAL N= %d E= %d DATAPACKING=BLOCK\n",g->nnodes/g->nstrand,g->nbfaces); 
  fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=NODAL, 4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED, 7=CELLCENTERED, 8=CELLCENTERED)\n");
  for(i=0;i<g->nnodes/g->nstrand;i++)
    fprintf(fp,"%f\n",g->x[3*i]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes/g->nstrand;i++)
    fprintf(fp,"%f\n",g->x[3*i+1]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes/g->nstrand;i++)
    fprintf(fp,"%f\n",g->x[3*i+2]);
  fprintf(fp,"\n");
  // rho
  for(i=0;i<g->nbfaces;i++)
    {
     rho = s->q[5*i];
     fprintf(fp,"%f\n",rho);
    }
  // u
  for(i=0;i<g->nbfaces;i++)
    {
     rho = s->q[5*i];
     rhou = s->q[5*i+1];
     fprintf(fp,"%f\n",rhou/rho);
    }
  // v
  for(i=0;i<g->nbfaces;i++)
    {
     rho = s->q[5*i];
     rhov = s->q[5*i+2];
     fprintf(fp,"%f\n",rhov/rho);
    }
  // w
  for(i=0;i<g->nbfaces;i++)
    {
     rho = s->q[5*i];
     rhow = s->q[5*i+3];
     fprintf(fp,"%f\n",rhow/rho);
    }
  // cp
  for(i=0;i<g->nbfaces;i++)
      {
       rho = s->q[5*i];
       rhou = s->q[5*i+1];
       rhov = s->q[5*i+2];
       rhow = s->q[5*i+3];
       e = s->q[5*i+4];
       pp=(gamm-1)*(e-0.5*(rhou*rhou+rhov*rhov+rhow*rhow)/rho);
       cp=(pp-pinf)/(0.5*s->mach*s->mach);  
       fprintf(fp,"%f\n",cp);
      }
  fprintf(fp,"\n");
  for(i=0;i<g->nbfaces;i++)
    fprintf(fp,"%d %d %d %d\n",g->conn[8*i],g->conn[8*i+1],g->conn[8*i+2],g->conn[8*i+3]);

  fclose(fp);
  }


  //3D wing
  if(g->test == 2)
  {
fprintf(fp,"VARIABLES = \"X\",\"Y\",\"Z\",\"RHO\",\"RHOU\",\"RHOV\",\"RHOW\",\"E\"\n");
fprintf(fp,"ZONE ZONETYPE=FEQUADRILATERAL N= %d E= %d DATAPACKING=BLOCK\n",g->nbfaces*4,g->nbfaces); 
fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=NODAL, 4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED, 7=CELLCENTERED, 8=CELLCENTERED)\n");
  
  for(i=0;i<g->nbfaces;i++) // x 
  {
     iface = g->bfaces[i];
     for(j=0;j<4;j++)
     {
       n = g->faces[8*iface+j];
       fprintf(fp,"%f\n",g->x[3*n]);
     }
  }
       fprintf(fp,"\n");
  for(i=0;i<g->nbfaces;i++)  // y
  {
     iface = g->bfaces[i];
     for(j=0;j<4;j++)
     {
       n = g->faces[8*iface+j];
       fprintf(fp,"%f\n",g->x[3*n+1]);
     }
  }
       fprintf(fp,"\n");
  for(i=0;i<g->nbfaces;i++) // z
  {
     iface = g->bfaces[i];
     for(j=0;j<4;j++)
     {
       n = g->faces[8*iface+j];
       fprintf(fp,"%f\n",g->x[3*n+2]);
     }
  }
       fprintf(fp,"\n");
  //flow variables
  for(n=0;n<NVAR;n++)
    for(i=0;i<g->nbfaces;i++)
      {
       iface = g->bfaces[i];
       icell = g->faces[8*iface+4];
       fprintf(fp,"%f\n",s->q[5*icell+n]);
      }
  fprintf(fp,"\n");

 for(i=0;i<g->nbfaces;i++) //connectivity
 {   
    fprintf(fp,"%d %d %d %d\n",4*i+1,4*i+2,4*i+3,4*i+4);
 }
  fclose(fp);
 }

// ==================================================================
//
// modifications for multiple processors (stitch for tecplot)
// 
// ==================================================================
  int     ii,k,totNumNode,totNumCell,totNumFace,stride,temp,temp2;
  int    *allNumNode    = (int *) malloc(sizeof(int)*nproc);
  int    *allNumNodeCum = (int *) malloc(sizeof(int)*nproc);
  int    *allNumCell    = (int *) malloc(sizeof(int)*nproc);
  int    *allNumCellCum = (int *) malloc(sizeof(int)*nproc);
  int    *allNumFace    = (int *) malloc(sizeof(int)*nproc);
  int    *allNumFaceCum = (int *) malloc(sizeof(int)*nproc);
  int    *myConn, *allConn, *surfConn, *tempConn;
  int    *subdomainconn, *subdomainOtherID;
  double *allPos, *allSoln;

  MPI_Allgather(&g->nnodes,  1, MPI_INT, allNumNode, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(&g->ncells,  1, MPI_INT, allNumCell, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(&g->nbfaces, 1, MPI_INT, allNumFace, 1, MPI_INT, MPI_COMM_WORLD);

  // compute the total number of nodes and cells (by a single proc)
  if (myid == 0)
  {
    totNumNode = totNumCell = totNumFace = 0;
    for (i = 0; i < nproc; i++)
    {
      totNumNode += allNumNode[i];
      totNumCell += allNumCell[i];
      totNumFace += allNumFace[i];
    }
  }
  

  // broadcast totNumNode and totNumCell to all procs (needed?)
  MPI_Bcast(&totNumNode,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&totNumCell,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&totNumFace,1,MPI_INT,0,MPI_COMM_WORLD);

  allNumNodeCum[0] = allNumCellCum[0] = allNumFaceCum[0] = 0;
  for (i = 1; i < nproc; i++)
  {
    allNumNodeCum[i] = allNumNodeCum[i-1] + allNumNode[i-1];
    allNumCellCum[i] = allNumCellCum[i-1] + allNumCell[i-1];
    allNumFaceCum[i] = allNumFaceCum[i-1] + allNumFace[i-1];
  }

  // Allgatherv the pos, soln, and conn values
  allPos = (double *) malloc(sizeof(double)*3*totNumNode);

  int *recvcounts   = (int *) malloc(sizeof(int)*nproc);
  int *displacement = (int *) malloc(sizeof(int)*nproc);

  //
  // allPos gather
  //
  for (i = 0; i < nproc; i++)
    recvcounts[i] = 3*allNumNode[i];

  displacement[0] = 0;
  for (i = 1; i < nproc; i++)
    displacement[i] = displacement[i-1] + recvcounts[i-1];

  MPI_Allgatherv(g->x, 3*g->nnodes, MPI_DOUBLE, allPos,
    recvcounts, displacement, MPI_DOUBLE, MPI_COMM_WORLD);

  //
  // allSoln gather
  //
  allSoln = (double *) malloc(sizeof(double)*5*totNumCell);
  for (i = 0; i < nproc; i++)
    recvcounts[i] = 5*allNumCell[i];

  displacement[0] = 0;
  for (i = 1; i < nproc; i++)
    displacement[i] = displacement[i-1] + recvcounts[i-1];

  MPI_Allgatherv(s->q, 5*g->ncells, MPI_DOUBLE, allSoln,
    recvcounts, displacement, MPI_DOUBLE, MPI_COMM_WORLD);

  //
  // allConn gather
  //
  myConn   = (int *) malloc(sizeof(int)*8*g->ncells);
  allConn  = (int *) malloc(sizeof(int)*8*totNumCell);
  surfConn = (int *) malloc(sizeof(int)*4*totNumFace);
  tempConn = (int *) malloc(sizeof(int)*4*totNumFace);

  temp = 8*g->ncells;
  for (i = 0; i < temp; i++)
  {
    myConn[i] = g->conn[i] + allNumNodeCum[myid];
  }

  for (i = 0; i < nproc; i++)
    recvcounts[i] = 8*allNumCell[i];

  displacement[0] = 0;
  for (i = 1; i < nproc; i++)
    displacement[i] = displacement[i-1] + recvcounts[i-1];

  MPI_Allgatherv(myConn, 8*g->ncells, MPI_INT, allConn,
    recvcounts, displacement, MPI_INT, MPI_COMM_WORLD);

  int *mySurfConn = (int *) malloc(sizeof(int)*4*g->nbfaces);
  int *tempMySurf = (int *) malloc(sizeof(int)*4*g->nbfaces);
  j = 0;
  for (i = 0; i < g->nbfaces; i++)
  {
    mySurfConn[j  ] = g->conn[8*i  ];
    mySurfConn[j+1] = g->conn[8*i+1];
    mySurfConn[j+2] = g->conn[8*i+2];
    mySurfConn[j+3] = g->conn[8*i+3];
    j+=4;
  }

  temp = 4*g->nbfaces;
  for (i = 0; i < temp; i++)
  {
    tempMySurf[i]  = mySurfConn[i] + allNumNodeCum[myid];
    mySurfConn[i] += allNumNodeCum[myid]/g->nstrand;
    
  }

  for (i = 0; i < nproc; i++)
    recvcounts[i] = 4*allNumFace[i];

  displacement[0] = 0;
  for (i = 1; i < nproc; i++)
    displacement[i] = displacement[i-1] + recvcounts[i-1];

  MPI_Allgatherv(mySurfConn, 4*g->nbfaces, MPI_INT, surfConn,
    recvcounts, displacement, MPI_INT, MPI_COMM_WORLD);
    
  MPI_Allgatherv(tempMySurf, 4*g->nbfaces, MPI_INT, tempConn,
    recvcounts, displacement, MPI_INT, MPI_COMM_WORLD);    

// ==================================================================
// Increase the conn number for subsequent domains and update
// the conn values based on subdomainconn.dat (from meshgen)
// ==================================================================
  if (myid == 0)
  {
    int nodeID;
    // rewrite connvalues based on subdomainconn.dat
    subdomainconn    = (int *) malloc(sizeof(int)*totNumNode);
    subdomainOtherID = (int *) malloc(sizeof(int)*totNumNode);

    fp = fopen("./QuadData/domain000/subdomainconn.dat","r");
    
    if (fp == NULL)
    {
      printf("WARNING: subdomainconn.dat file missing. Recheck. Stopping.\n");
      MPI_Abort(MPI_COMM_WORLD,33);
    }
    fscanf(fp,"%*d"); // ignore the number
    for (i = 0; i < totNumNode; i++)
      fscanf(fp,"%d",&subdomainconn[i]);

    fscanf(fp,"%*d"); // ignore the number
    for (i = 0; i < totNumNode; i++)
    {
      fscanf(fp,"%d",&subdomainOtherID[i]);
    }

    fclose(fp);

    //
    // allconn 
    //
    for (i = 0; i < totNumCell; i++)
    {
      for (j = 0; j < 8; j++)
      {
        temp   = 8*i+j;
        nodeID = allConn[temp]-1;

        if (subdomainconn[nodeID] != -1)
          allConn[temp] = subdomainconn[nodeID]+1;

      }
    }
    //
    // surfconn
    //
    // adhoc fix variables
    int flag,diff,id,domainid;
    double flr;

    k = 0;
    for (i = 0; i < totNumFace; i++)
    {
      for (j = 0; j < 4; j++)
      {
        id = subdomainconn[tempConn[k]-1];
        if (id != -1)
        {

          domainid = subdomainOtherID[tempConn[k]-1];

          flag = 0;
          if(domainid > 0)
          {
            for (ii = nproc-1; ii >= domainid; ii--)
            {
              flr = (double)(id/allNumNodeCum[ii]);
              if(flr == 1.0)
              {
                diff = id - allNumNodeCum[ii];
                flag = 1;
                break;
              }
            }
          }
          
          if(flag==0)
          {
            diff = id;
          }
          surfConn[k] = diff+1+allNumNodeCum[domainid]/g->nstrand;
        }
        k++;
      }
    }

// ==================================================================
// Write stitched outputs into file
// ==================================================================

    //
    // VOLUME FILE
    //
    fp=fopen("./output/volall.dat","w");

    fprintf(fp,"VARIABLES = \"X\",\"Y\",\"Z\",\"RHO\",\"U\",\"V\",\"W\",\"CP\"\n");
    fprintf(fp,"ZONE ZONETYPE=FEBRICK N= %d E= %d DATAPACKING=BLOCK\n",totNumNode,totNumCell); 
    fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=NODAL, 4=CELLCENTERED, 5=CELLCENTERED, "
        "6=CELLCENTERED, 7=CELLCENTERED, 8=CELLCENTERED)\n");

     // list of x-coordinate positions 
     for(i = 0; i < totNumNode; i++)
       fprintf(fp,"%f\n",allPos[3*i]);
     fprintf(fp,"\n");

     // list of y-coordinate positions
     for(i = 0; i < totNumNode; i++)
       fprintf(fp,"%f\n",allPos[3*i+1]);
     fprintf(fp,"\n");

     // list of a-coordinate positions
     for(i = 0; i < totNumNode; i++)
       fprintf(fp,"%f\n",allPos[3*i+2]);
     fprintf(fp,"\n");

    //rho
    for(i = 0; i < totNumCell; i++)
    {
      rho = allSoln[5*i];
      fprintf(fp,"%f\n",rho);
    }
    // u
    for(i=0; i < totNumCell; i++)
    {
      rho  = allSoln[5*i];
      rhou = allSoln[5*i+1];
      fprintf(fp,"%f\n",rhou/rho);
    }
    // v
    for(i=0; i < totNumCell; i++)
    {
      rho  = allSoln[5*i];
      rhov = allSoln[5*i+2];
      fprintf(fp,"%f\n",rhov/rho);
    }
    // w
    for(i=0; i < totNumCell; i++)
    {
      rho  = allSoln[5*i];
      rhow = allSoln[5*i+3];
      fprintf(fp,"%f\n",rhow/rho);
    }
    // cp
    for(i=0; i < totNumCell; i++)
    {
      rho  = allSoln[5*i];
      rhou = allSoln[5*i+1];
      rhov = allSoln[5*i+2];
      rhow = allSoln[5*i+3];
      e    = allSoln[5*i+4];
      pp   = (gamm-1)*(e-0.5*(rhou*rhou+rhov*rhov+rhow*rhow)/rho);
      cp   = (pp-pinf)/(0.5*s->mach*s->mach); 
      fprintf(fp,"%f\n",cp);
    }
    fprintf(fp,"\n");

    // connectivity information for the cells
    for(i = 0; i < totNumCell; i++)
    {
      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
      allConn[8*i  ],allConn[8*i+1],allConn[8*i+2],allConn[8*i+3],
      allConn[8*i+4],allConn[8*i+5],allConn[8*i+6],allConn[8*i+7]);
    }

    fclose(fp);

    // 
    // SURFACE FILE
    //
    temp = totNumNode/g->nstrand;

    int *surfid = (int *) malloc(sizeof(int)*temp);

    k = 0;
    for (i = 0 ; i < nproc; i++)
    {
      temp2 = allNumNode[i]/g->nstrand;
      for (j = 0; j < temp2; j++)
      {
        surfid[k] = allNumNodeCum[i] + j;
        k++;
      } // j loop
    } // i loop

    fp=fopen("./output/surfall.dat","w");

    fprintf(fp,"VARIABLES = \"X\",\"Y\",\"Z\",\"RHO\",\"U\",\"V\",\"W\",\"CP\"\n");
    fprintf(fp,"ZONE ZONETYPE=FEQUADRILATERAL N= %d E= %d DATAPACKING=BLOCK\n",temp,totNumFace); 
    fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=NODAL, 4=CELLCENTERED, 5=CELLCENTERED, "
        "6=CELLCENTERED, 7=CELLCENTERED, 8=CELLCENTERED)\n");
    // x
    for(i = 0; i < temp; i++)
      fprintf(fp,"%f\n",allPos[3*surfid[i]  ]);
    fprintf(fp,"\n");
    // y
    for(i = 0; i < temp; i++)
      fprintf(fp,"%f\n",allPos[3*surfid[i]+1]);
    fprintf(fp,"\n");
    // z
    for(i = 0; i < temp; i++)
      fprintf(fp,"%f\n",allPos[3*surfid[i]+2]);
    fprintf(fp,"\n");

    k = 0;
    for (i = 0 ; i < nproc; i++)
    {
      temp2 = allNumCell[i]/(g->nstrand-1);
      for (j = 0; j < temp2; j++)
      {
        surfid[k] = allNumCellCum[i] + j;
        k++;
      } // j loop
    } // i loop

    // rho
    for(i = 0; i < totNumFace; i++)
      {
       rho = allSoln[5*surfid[i]];
       fprintf(fp,"%f\n",rho);
      }
    // u
    for(i = 0; i < totNumFace; i++)
      {
       rho  = allSoln[5*surfid[i]];
       rhou = allSoln[5*surfid[i]+1];
       fprintf(fp,"%f\n",rhou/rho);
      }
    // v
    for(i = 0; i < totNumFace; i++)
      {
       rho  = allSoln[5*surfid[i]];
       rhov = allSoln[5*surfid[i]+2];
       fprintf(fp,"%f\n",rhov/rho);
      }
    // w
    for(i = 0; i < totNumFace; i++)
      {
       rho  = allSoln[5*surfid[i]];
       rhow = allSoln[5*surfid[i]+3];
       fprintf(fp,"%f\n",rhow/rho);
      }
    // cp
    for(i = 0; i < totNumFace; i++)
        {
         rho  = allSoln[5*surfid[i]  ];
         rhou = allSoln[5*surfid[i]+1];
         rhov = allSoln[5*surfid[i]+2];
         rhow = allSoln[5*surfid[i]+3];
         e    = allSoln[5*surfid[i]+4];
         pp=(gamm-1)*(e-0.5*(rhou*rhou+rhov*rhov+rhow*rhow)/rho);
         cp=(pp-pinf)/(0.5*s->mach*s->mach);  
         fprintf(fp,"%f\n",cp);
        }
    fprintf(fp,"\n");


    for(i = 0; i  < totNumFace; i++)
    {
      fprintf(fp,"%d %d %d %d\n",
        surfConn[4*i],surfConn[4*i+1],surfConn[4*i+2],surfConn[4*i+3]);
    }

    fclose(fp);

  }

}



// ##################################################################
//
// outputdq
//
// ##################################################################
void outputdq(GRID *g,SOLN *s)
{
  int i,n;
  int iface,node1,node2;
  int icell;
  double x1,y1,x2,y2,rho,rhou,rhov,e,pp,cp;
  FILE *fp;
  char fname[80];
  static int istep0=0;
  //
  sprintf(fname,"dq%d.plt",istep0);
  fp=fopen(fname,"w");
  fprintf(fp,"VARIABLES = \"X\",\"Y\",\"RHO\",\"RHOU\",\"RHOV\",\"E\"\n");
  fprintf(fp,"ZONE ZONETYPE=FEQUADRILATERAL N=%d E=%d DATAPACKING=BLOCK\n",g->nnodes,g->ncells); 
  fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=CELLCENTERED, 4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED)\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[2*i]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[2*i+1]);
  fprintf(fp,"\n");
  for(n=0;n<NVAR;n++)
    for(i=0;i<g->ncells;i++)
      fprintf(fp,"%f\n",s->ddq[4*i+n]);
  fprintf(fp,"\n");
  for(i=0;i<g->ncells;i++)
    fprintf(fp,"%d %d %d %d\n",g->conn[4*i]+1,g->conn[4*i+1]+1,g->conn[4*i+2]+1,g->conn[4*i+3]+1);
  fclose(fp);
  istep0++;
}

void outputr(GRID *g,SOLN *s)
{
  int i,n;
  int iface,node1,node2;
  int icell;
  double x1,y1,x2,y2,rho,rhou,rhov,e,pp,cp;
  FILE *fp;
  char fname[80];
  static int istep0=0;
  //
  sprintf(fname,"r%d.plt",istep0);
  fp=fopen(fname,"w");
  fprintf(fp,"VARIABLES = \"X\",\"Y\",\"RHO\",\"RHOU\",\"RHOV\",\"E\"\n");
  fprintf(fp,"ZONE ZONETYPE=FEQUADRILATERAL N=%d E=%d DATAPACKING=BLOCK\n",g->nnodes,g->ncells); 
  fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=CELLCENTERED, 4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED)\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[2*i]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[2*i+1]);
  fprintf(fp,"\n");
  for(n=0;n<NVAR;n++)
    for(i=0;i<g->ncells;i++)
      fprintf(fp,"%f\n",s->r[4*i+n]);
  fprintf(fp,"\n");
  for(i=0;i<g->ncells;i++)
    fprintf(fp,"%d %d %d %d\n",g->conn[4*i]+1,g->conn[4*i+1]+1,g->conn[4*i+2]+1,g->conn[4*i+3]+1);
  fclose(fp);
  istep0++;
}


