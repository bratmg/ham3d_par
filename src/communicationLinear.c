  #include <stdio.h>
  #include "ham2dtypes.h"
  #include "ham2dFunctionDefs.h"
  #include <stdlib.h>
  #define NQ 5
  #include <mpi.h>

  void communicationLinear(GRID *g, SOLN *s,int myid)
  {

  int i,j,iq,ic;
  int iadjp,ist,ilengt,iend;
  int tag,error,ierr;
  int c1,c2,c3;
  int order;

  double *dq1r,*dq2r,*dq3r,*dq4r,*dq5r;
  double *dq1s,*dq2s,*dq3s,*dq4s,*dq5s;
  double *dq1ra,*dq2ra,*dq3ra,*dq4ra,*dq5ra;
  double *dq1sa,*dq2sa,*dq3sa,*dq4sa,*dq5sa;
  double *dq1rb,*dq2rb,*dq3rb,*dq4rb,*dq5rb;
  double *dq1sb,*dq2sb,*dq3sb,*dq4sb,*dq5sb;

  order = g->order;
  MPI_Request request[15];
  MPI_Status status[15];

  dq1r = (double *) malloc(sizeof(double)*g->ncommu);
  dq2r = (double *) malloc(sizeof(double)*g->ncommu);
  dq3r = (double *) malloc(sizeof(double)*g->ncommu);
  dq4r = (double *) malloc(sizeof(double)*g->ncommu);
  dq5r = (double *) malloc(sizeof(double)*g->ncommu);


  dq1s = (double *) malloc(sizeof(double)*g->ncommu);
  dq2s = (double *) malloc(sizeof(double)*g->ncommu);
  dq3s = (double *) malloc(sizeof(double)*g->ncommu);
  dq4s = (double *) malloc(sizeof(double)*g->ncommu);
  dq5s = (double *) malloc(sizeof(double)*g->ncommu);
 
  dq1ra = (double *) malloc(sizeof(double)*g->ncommu);
  dq2ra = (double *) malloc(sizeof(double)*g->ncommu);
  dq3ra = (double *) malloc(sizeof(double)*g->ncommu);
  dq4ra = (double *) malloc(sizeof(double)*g->ncommu);
  dq5ra = (double *) malloc(sizeof(double)*g->ncommu);


  dq1sa = (double *) malloc(sizeof(double)*g->ncommu);
  dq2sa = (double *) malloc(sizeof(double)*g->ncommu);
  dq3sa = (double *) malloc(sizeof(double)*g->ncommu);
  dq4sa = (double *) malloc(sizeof(double)*g->ncommu);
  dq5sa = (double *) malloc(sizeof(double)*g->ncommu);

  if(order==5)
  {
  dq1rb = (double *) malloc(sizeof(double)*g->ncommu);
  dq2rb = (double *) malloc(sizeof(double)*g->ncommu);
  dq3rb = (double *) malloc(sizeof(double)*g->ncommu);
  dq4rb = (double *) malloc(sizeof(double)*g->ncommu);
  dq5rb = (double *) malloc(sizeof(double)*g->ncommu);


  dq1sb = (double *) malloc(sizeof(double)*g->ncommu);
  dq2sb = (double *) malloc(sizeof(double)*g->ncommu);
  dq3sb = (double *) malloc(sizeof(double)*g->ncommu);
  dq4sb = (double *) malloc(sizeof(double)*g->ncommu);
  dq5sb = (double *) malloc(sizeof(double)*g->ncommu);
  }


  //copy to psil
  for(i=0;i<g->nadjp;i++)
  {
    iadjp  = g->iadjp[i];
    ist    = g->istp[i];
    ilengt = g->ilengp[i];
    iend   = ist + ilengt;
    
    for(j=ist;j<iend;j++)
    { 
 
       c2 = g->irecvconn2[j*6+1];
       c1 = g->irecvconn2[j*6+2];

      
       
       if(order==5)
       {
       c3 = g->irecvconn2[j*6]; 
       
       // third cell
       dq1sb[j] = s->dq[NVAR*c1];
       dq2sb[j] = s->dq[NVAR*c1+1];
       dq3sb[j] = s->dq[NVAR*c1+2];
       dq4sb[j] = s->dq[NVAR*c1+3];
       dq5sb[j] = s->dq[NVAR*c1+4];
       }


       // second cell 
       dq1s[j] = s->dq[NVAR*c2];
       dq2s[j] = s->dq[NVAR*c2+1];
       dq3s[j] = s->dq[NVAR*c2+2];
       dq4s[j] = s->dq[NVAR*c2+3];
       dq5s[j] = s->dq[NVAR*c2+4];

       // first cell
       dq1sa[j] = s->dq[NVAR*c1];
       dq2sa[j] = s->dq[NVAR*c1+1];
       dq3sa[j] = s->dq[NVAR*c1+2];
       dq4sa[j] = s->dq[NVAR*c1+3];
       dq5sa[j] = s->dq[NVAR*c1+4];

    } 
  }


  for(i=0;i<g->nadjp;i++)
  {
  iadjp  = g->iadjp[i];
  ist    = g->istp[i];
  ilengt = g->ilengp[i];

  tag = 15*(2*(myid+1)+2*(iadjp+1));

  error = MPI_Irecv(&dq1r[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[0]);
  tag = tag+1;
  error = MPI_Irecv(&dq2r[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[1]);
tag = tag+1;
  error = MPI_Irecv(&dq3r[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[2]);
tag = tag+1;
  error = MPI_Irecv(&dq4r[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[3]);
  
tag = tag+1;
  error = MPI_Irecv(&dq5r[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[4]);
  
tag = tag+1;
   error = MPI_Irecv(&dq1ra[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[5]);
  tag = tag +1;
  error = MPI_Irecv(&dq2ra[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[6]);
tag = tag +1;
  error = MPI_Irecv(&dq3ra[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[7]);
tag = tag +1;
  error = MPI_Irecv(&dq4ra[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[8]);
tag = tag +1;
  error = MPI_Irecv(&dq5ra[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[9]);

if(order==5) {
tag = tag+1;
  error = MPI_Irecv(&dq1rb[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[10]);
  tag = tag +1;
  error = MPI_Irecv(&dq2rb[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[11]);
tag = tag +1;
  error = MPI_Irecv(&dq3rb[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[12]);
tag = tag +1;
  error = MPI_Irecv(&dq4rb[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[13]);
tag = tag +1;
  error = MPI_Irecv(&dq5rb[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[14]);
}
}
  
  
  for(i=0;i<g->nadjp;i++)
  {
  iadjp  = g->iadjp[i];
  ist    = g->istp[i];
  ilengt = g->ilengp[i];
  tag = 15*(2*(myid+1)+2*(iadjp+1));
 

  error = MPI_Send(&dq1s[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD);
  tag = tag+1;
  error = MPI_Send(&dq2s[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD);
tag = tag+1;
  error = MPI_Send(&dq3s[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD);
tag = tag+1;
  error = MPI_Send(&dq4s[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD);
 tag = tag+1;
  error = MPI_Send(&dq5s[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD);
 
 tag = tag+1;
  error = MPI_Send(&dq1sa[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
tag = tag +1;

  error = MPI_Send(&dq2sa[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
tag = tag +1;

  error = MPI_Send(&dq3sa[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
tag = tag +1;

error = MPI_Send(&dq4sa[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
tag = tag +1;
error = MPI_Send(&dq5sa[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);

if(order==5){
tag = tag +1;

  error = MPI_Send(&dq1sb[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
tag = tag +1;

  error = MPI_Send(&dq2sb[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
tag = tag +1;

  error = MPI_Send(&dq3sb[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
tag = tag +1;

  error = MPI_Send(&dq4sb[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
  tag = tag +1;
  error = MPI_Send(&dq5sb[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
}
  }
  
  for(i=0;i<g->nadjp;i++)
  {
    error=MPI_Wait(&request[0],&status[0]);
    error=MPI_Wait(&request[1],&status[1]);
    error=MPI_Wait(&request[2],&status[2]);
    error=MPI_Wait(&request[3],&status[3]);
    error=MPI_Wait(&request[4],&status[4]);
    error=MPI_Wait(&request[5],&status[5]);
    error=MPI_Wait(&request[6],&status[6]);
    error=MPI_Wait(&request[7],&status[7]);
    error=MPI_Wait(&request[8],&status[8]);
    error=MPI_Wait(&request[9],&status[9]);
    if(order==5){
    error=MPI_Wait(&request[10],&status[10]);
    error=MPI_Wait(&request[11],&status[11]);
    error=MPI_Wait(&request[12],&status[12]);
    error=MPI_Wait(&request[13],&status[13]);
    error=MPI_Wait(&request[14],&status[14]);
    }
  }

 ierr = MPI_Barrier(MPI_COMM_WORLD);
 
  // copy from receive data to designated cell 
  for(i=0;i<g->nadjp;i++)
  {
    iadjp  = g->iadjp[i];
    ist    = g->istp[i];
    ilengt = g->ilengp[i];
    iend   = ist + ilengt;
    for(j=ist;j<iend;j++)
    {
      c2 = g->irecv[j*2]; // second cell index            
      c1 = g->irecvconn[j];
      
      // for duplication case
      if(g->idup[c2]==2)
      {
        g->dpsil_deri[c1][0]   = dq1r[j]; 
        g->dpsil_deri[c1][1]   = dq2r[j];   
        g->dpsil_deri[c1][2]   = dq3r[j];   
        g->dpsil_deri[c1][3]   = dq4r[j];
        g->dpsil_deri[c1][4]   = dq5r[j];

        g->dpsila_deri[c1][0]  = dq1ra[j]; 
        g->dpsila_deri[c1][1]  = dq2ra[j];   
        g->dpsila_deri[c1][2]  = dq3ra[j];   
        g->dpsila_deri[c1][3]  = dq4ra[j];
        g->dpsila_deri[c1][4]  = dq5ra[j];

        if(order==5){
        g->dpsilb_deri[c1][0]  = dq1rb[j]; 
        g->dpsilb_deri[c1][1]  = dq2rb[j];   
        g->dpsilb_deri[c1][2]  = dq3rb[j];   
        g->dpsilb_deri[c1][3]  = dq4rb[j];
        g->dpsilb_deri[c1][4]  = dq5rb[j];
        }
      }
      else
      {
        g->psil_deri[c2][0]  = dq1r[j]; 
        g->psil_deri[c2][1]  = dq2r[j];   
        g->psil_deri[c2][2]  = dq3r[j];   
        g->psil_deri[c2][3]  = dq4r[j];
        g->psil_deri[c2][4]  = dq4r[j];

        g->psila_deri[c2][0]  = dq1ra[j]; 
        g->psila_deri[c2][1]  = dq2ra[j];   
        g->psila_deri[c2][2]  = dq3ra[j];   
        g->psila_deri[c2][3]  = dq4ra[j];
        g->psila_deri[c2][4]  = dq4ra[j];

        if(order==5){
        g->psilb_deri[c2][0]  = dq1rb[j]; 
        g->psilb_deri[c2][1]  = dq2rb[j];   
        g->psilb_deri[c2][2]  = dq3rb[j];   
        g->psilb_deri[c2][3]  = dq4rb[j];
        g->psilb_deri[c2][4]  = dq4rb[j];
        }
      }

  //if(myid==0)
  //{
  //printf("myid:%d, q1r[j]:%f, q1s[j]:%f\n",myid,dq1r[j],dq1s[j]);
  //printf("myid:%d, q2r[j]:%f, q2s[j]:%f\n",myid,dq2r[j],dq2s[j]);
  //printf("myid:%d, q3r[j]:%f, q3s[j]:%f\n",myid,dq3r[j],dq3s[j]);
  //printf("myid:%d, q4r[j]:%f, q4s[j]:%f\n",myid,dq4r[j],dq4s[j]);
  //}

    } 
  }
 

  //if(myid==1)
  //{
  //printf("myid:%d, q1r[5]:%f, q1s[5]:%f\n",myid,q1r[0],q1s[0]);
  //printf("myid:%d, q1r[7]:%f, q1s[7]:%f\n",myid,q1r[1],q1s[1]);
  //printf("myid:%d, q1r[9]:%f, q1s[9]:%f\n",myid,q1r[2],q1s[2]);
  //printf("myid:%d, q1r[11]:%f, q1s[11]:%f\n",myid,q1r[3],q1s[3]);
  //}

  free(dq1r);
  free(dq2r);
  free(dq3r);
  free(dq4r);
  free(dq5r);

  free(dq1s);
  free(dq2s);
  free(dq3s);
  free(dq4s);
  free(dq5s);
  
  free(dq1ra);
  free(dq2ra);
  free(dq3ra);
  free(dq4ra);
  free(dq5ra);

  free(dq1sa);
  free(dq2sa);
  free(dq3sa);
  free(dq4sa);
  free(dq5sa);

  if(order==5)
  {
  free(dq1rb);
  free(dq2rb);
  free(dq3rb);
  free(dq4rb);
  free(dq5rb);

  free(dq1sb);
  free(dq2sb);
  free(dq3sb);
  free(dq4sb);
  free(dq5sb);
  }


  // subroutine COMMUNICATION end
  }
