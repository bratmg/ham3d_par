// HAM3D V6
// 1. strand grid is set
// 2. Inviscid solid surface boundary condition is set
// 3. Only explicit time integration can be used
// 4. It does not need nearbody Radius
// 5. periodic boundary condition is set
// 6. Add output option
// 7. integrate 3D wing case 2014.12.17

//======================================================
//bug fix report
//
//1. computeRHS.c : line 357 (2014.11.19)
// else if (rightCell==-2) 
// -> else if(rightCell==-2 && g->test==0)
//2. computeRHS.c: line339
// if(rightCell>-1||g->test==1)->if(rightCell>-1)
//
//3. same as computeRHSk.c : line 310, line:325, line: 396
//4. outputSolution.c fix
//4. apply_periodic_LHS.c bug fixed : 2014.11.21
//5. jac_roe.f90 : add the muliply at last line
//6. error fixed in apply_periodic_LHS.c : 2014.12.17
//======================================================

#include <stdio.h>
#include <stdlib.h>
#include "ham3dtypes.h"
#include "ham3dFunctionDefs.h"
#include <time.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char **argv)
{
  GRID *g;
  SOLN *s;
  int ngrids;
  int nsteps;
  int nwrite;
  int i,n,nn;
  double dt;
  double CFL;
  double l2rho,linfrho;
  double sref;
  int msweep;
  char c;
  char fname[80];
  char scheme[20];
  char timeInteg[20];
  int order,timeacc,visc,test,outform;
  FILE *fp;
  clock_t start, end;
  double cpu_time_used;
  // 
  int nproc,myid,ierr;

  ierr = MPI_Init(&argc,&argc);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  printf("Hello world! I'm process %i out of %i process\n",myid,nproc);
  ierr = MPI_Barrier(MPI_COMM_WORLD);

  if(myid==0)welcome();

  // Read in file input from input.ham2d
  fp=fopen("input.ham3d","r");
  while((c=fgetc(fp))!='\n');
  while((c=fgetc(fp))!='\n');
  while((c=fgetc(fp))!='\n');
  while((c=fgetc(fp))!='\n');
  fscanf(fp,"scheme=%s\n",scheme);
  fscanf(fp,"time integration=%s\n",timeInteg);
  fscanf(fp,"order=%d\n",&order);
  fscanf(fp,"timeacc=%d\n",&timeacc);
  fscanf(fp,"nsteps=%d\n",&nsteps);
  fscanf(fp,"nwrite=%d\n",&nwrite);
  fscanf(fp,"dt=%lf\n",&dt);
  fscanf(fp,"CFL=%lf\n",&CFL);
  fscanf(fp,"msweep=%d\n",&msweep);
  fscanf(fp,"visc=%d\n",&visc);
  fscanf(fp,"testcase=%d\n",&test);
  fscanf(fp,"output=%d\n",&outform);
  fclose(fp);

  // if(myid==0)
  // {
  // trace(nsteps);
  // tracef(dt);
  // }

  //
  ngrids=1;
  g=(GRID *) malloc(sizeof(GRID)*ngrids);
  s=(SOLN *) malloc(sizeof(SOLN)*ngrids);
  //
  // preprocess grids
  // code is written with an overset
  // framework in mind (but not implemented yet)
  // (ngrid==1) for now
  //

  for(i=0;i<ngrids;i++) 
    {
      readGrid(&g[i],myid,nproc); //need to strand grid reading and make connectivity
      g[i].visc=visc;
      g[i].test=test;
      s[i].scheme = scheme;
      s[i].outform=outform;
      s[i].nsteps = nsteps;
      
      preprocess(&g[i],myid);
      initflow(&g[i],&s[i],myid);

      g[i].order=order;
      g[i].timeacc=timeacc;
      g[i].CFL=CFL;  
      s[i].cflnum=CFL;
     // s[i].cflmax=CFL;

      g[i].msweep=msweep;

      if (strcmp(timeInteg,"bdf1")==0) g[i].timeInteg=BDF1;
      if (strcmp(timeInteg,"bdf2")==0) g[i].timeInteg=BDF2;
    }

  //nn = 1;
  //outputSolution(&g[0],&s[0],nn,myid);
  //ierr = MPI_Barrier(MPI_COMM_WORLD);
  //exit(1);

  if(myid==0) basicScreenOutput(g,s);
 
  // =========================================================
  // now run chosen number of time steps
  // =========================================================
  if(myid<10) 
  sprintf(fname,"./output/sol_his00%i.dat",myid);
  else if(myid<100)
  sprintf(fname,"./output/sol_his0%i.dat",myid);
  

  fp = fopen(fname,"w");

  // if(myid==0)
  // {
  //   printf("#ham3d : using %s scheme for inversion\n",scheme);  
  // }
  cpu_time_used=0;

  ierr = MPI_Barrier(MPI_COMM_WORLD);
  // main iteration
  for (n=0;n<nsteps;n++) 
  { 
    for(i=0;i<ngrids;i++) 
    { 
      start=clock();

      stepSolution(scheme,&g[i],&s[i],dt,&l2rho,&linfrho,myid);
      
      computeForce(&g[i],&s[i]);
      
      end = clock();
      
      cpu_time_used+=(((double) (end - start)) / CLOCKS_PER_SEC);
       
      //
      // print to screen
      //
      if(myid == 0 && n == 0)
      {
        printf("\n --------------------------------------------------------------------------------------------------\n");
        printf(" |   Iter   | L_2 (density) | L_inf (density) |    Cl (lift)   |     Cd (drag)   | CPU Time (sec) |\n");
        printf(" --------------------------------------------------------------------------------------------------\n");
      }

      if(myid == 0) printf(" |  %*d  |  %e |   %e  |  %e |   %e  |  %e  |\n",
         6,n,l2rho,linfrho,s[i].cl,s[i].cd,cpu_time_used );

      //if(myid==0) printf("%d %e %e %2.4f %2.4f %2.4f\n",n,l2rho,linfrho,s[i].cl,s[i].cd,cpu_time_used);
      fprintf(fp,"%d %e %e %2.4f %2.4f %2.4f\n",n,l2rho,linfrho,s[i].cl,s[i].cd,cpu_time_used);


      if(myid == 0 && n == nsteps)
      {
        printf("\n --------------------------------------------------------------------------------------------------\n");
        printf(" |   Iter   | L_2 (density) | L_inf (density) |    Cl (lift)   |     Cd (drag)   | CPU Time (sec) |\n");
        printf(" --------------------------------------------------------------------------------------------------\n");
      }

    }
  
    ierr = MPI_Barrier(MPI_COMM_WORLD);
        
    //write solution
    if((n+1)%nwrite==0||n==0)
    { 
      nn = (n+1)/nwrite;
      outputSolution(&g[0],&s[0],nn,myid,nproc); 
    }

  } // n loop
  fclose(fp);
  
  ierr = MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}

//######################################################################
// END OF PROGRAM
//######################################################################

