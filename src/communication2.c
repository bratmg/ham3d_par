  #include <stdio.h>
  #include "ham3dtypes.h"
  #include "ham3dFunctionDefs.h"
  #include <stdlib.h>
  #define NQ 5
  #include <mpi.h>

  void communication2(GRID *g, SOLN *s,int myid)
  {

  int i,j,iq,ic;
  int iadjp,ist,ilengt,iend;
  int tag,error,ierr;
  int c1,c2;
 
  double *q1ra,*q2ra,*q3ra,*q4ra,*q5ra;
  double *q1sa,*q2sa,*q3sa,*q4sa,*q5sa;

  MPI_Request request[10];
  MPI_Status status[10];

  q1ra = (double *) malloc(sizeof(double)*g->ncommu);
  q2ra = (double *) malloc(sizeof(double)*g->ncommu);
  q3ra = (double *) malloc(sizeof(double)*g->ncommu);
  q4ra = (double *) malloc(sizeof(double)*g->ncommu);
  q5ra = (double *) malloc(sizeof(double)*g->ncommu);

  q1sa = (double *) malloc(sizeof(double)*g->ncommu);
  q2sa = (double *) malloc(sizeof(double)*g->ncommu);
  q3sa = (double *) malloc(sizeof(double)*g->ncommu);
  q4sa = (double *) malloc(sizeof(double)*g->ncommu);
  q5sa = (double *) malloc(sizeof(double)*g->ncommu);

  //copy to psil
  for(i=0;i<g->nadjp;i++)
  {
    iadjp  = g->iadjp[i];
    ist    = g->istp[i];
    ilengt = g->ilengp[i];
    iend   = ist + ilengt;
    
    for(j=ist;j<iend;j++)
    { 
       iq = g->isend[j*2]; // second neighbor cell index     
       c2 = g->qconni[iq*2];  // second cell index 
       c1 = g->qconni[iq*2+1];  // second cell index 
       
       //first cell
       q1sa[j] = s->q[NVAR*c1];
       q2sa[j] = s->q[NVAR*c1+1];
       q3sa[j] = s->q[NVAR*c1+2];
       q4sa[j] = s->q[NVAR*c1+3];
       q5sa[j] = s->q[NVAR*c1+4];

    } 
  }


  for(i=0;i<g->nadjp;i++)
  {
  iadjp  = g->iadjp[i];
  ist    = g->istp[i];
  ilengt = g->ilengp[i];

  tag = 10*(2*(myid+1)+2*(iadjp+1));

  error = MPI_Irecv(&q1ra[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[5]);
  tag = tag +1;
  
  error = MPI_Irecv(&q2ra[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[6]);
  tag = tag +1;

  error = MPI_Irecv(&q3ra[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[7]);
  tag = tag +1;

  error = MPI_Irecv(&q4ra[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[8]);
  tag = tag +1;

  error = MPI_Irecv(&q5ra[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[9]);

  }

  for(i=0;i<g->nadjp;i++)
  {
  iadjp  = g->iadjp[i];
  ist    = g->istp[i];
  ilengt = g->ilengp[i];
  
  tag = 10*(2*(myid+1)+2*(iadjp+1));

  error = MPI_Send(&q1sa[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
  tag = tag +1;

  error = MPI_Send(&q2sa[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
  tag = tag +1;

  error = MPI_Send(&q3sa[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
  tag = tag +1;

  error = MPI_Send(&q4sa[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
  tag = tag +1;

  error = MPI_Send(&q5sa[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);


  }
  
  for(i=0;i<g->nadjp;i++)
  {
    MPI_Wait(&request[5],&status[5]);
    MPI_Wait(&request[6],&status[6]);
    MPI_Wait(&request[7],&status[7]);
    MPI_Wait(&request[8],&status[8]);
    MPI_Wait(&request[9],&status[9]);

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
      c1 = g->irecvconn[j]; // first cell index
      
      if(g->idup[c2]==2)
      {
        g->dpsila[c1][0]  = q1ra[j]; 
        g->dpsila[c1][1]  = q2ra[j];   
        g->dpsila[c1][2]  = q3ra[j];   
        g->dpsila[c1][3]  = q4ra[j];
        g->dpsila[c1][4]  = q5ra[j];


      }
      else
      {
      g->psila[c2][0]  = q1ra[j]; 
      g->psila[c2][1]  = q2ra[j];   
      g->psila[c2][2]  = q3ra[j];   
      g->psila[c2][3]  = q4ra[j];
      g->psila[c2][4]  = q5ra[j];


      }
  //if(myid==0)
  //{
  //printf("myid:%d, q1r[j]:%f, q1s[j]:%f\n",myid,q1r[j],q1s[j]);
  //printf("myid:%d, q2r[j]:%f, q2s[j]:%f\n",myid,q2r[j],q2s[j]);
  //printf("myid:%d, q3r[j]:%f, q3s[j]:%f\n",myid,q3r[j],q3s[j]);
  //printf("myid:%d, q4r[j]:%f, q4s[j]:%f\n",myid,q4r[j],q4s[j]);
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


  free(q1ra);
  free(q2ra);
  free(q3ra);
  free(q4ra);
  free(q5ra);


  free(q1sa);
  free(q2sa);
  free(q3sa);
  free(q4sa);
  free(q5sa);


  
  // subroutine COMMUNICATION end
  }
