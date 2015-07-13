#include "ham3dtypes.h"
#include "ham3dFunctionDefs.h"

// ##################################################################
//
// WELCOME
//
// ##################################################################
void welcome(void)
{
   printf("\n");
   printf("######################################################################\n");
   printf("#                                                                    #\n");
   printf("#                            HAMSTR                                  #\n");
   printf("#                                                                    #\n");
   printf("#              HAMILTONIAN-STRAND FINITE VOLUME SOLVER               #\n");
   printf("#                                                                    #\n");
   printf("######################################################################\n");
}


// ##################################################################
//
// THANKS
//
// ##################################################################
void thanks(void)
{
   printf("\n");
   printf("######################################################################\n");
   printf("#                                                                    #\n");
   printf("#                       END OF SIMULATION                            #\n");
   printf("#                                                                    #\n");
   printf("######################################################################\n");
}

// ##################################################################
//
// BASICSCREENOUTPUT
//
// ##################################################################
void basicScreenOutput(GRID *g, SOLN *s)
{
   printf("======================================================================\n");
   printf("                       GRID RELATED OUTPUTS                           \n");
   printf("======================================================================\n");
   printf("\n");
   printf(" Number of nodes            %d\n",g->nnodes);
   printf(" Number of cells            %d\n",g->ncells);
   printf(" Number of faces            %d\n",g->nfaces);
   printf(" Number of colours          %d\n",g->ncolors);
   printf(" Number of chains           %d\n",g->nchains);
   printf(" Number of chain faces      %d\n",g->nchainFaces);
   printf(" Maximum length of chain    %d\n",g->nmaxchain);
   // printf(" Minimum cell volume        %0.8e\n",volmin);
   // printf(" Maximum cell volume        %0.8e\n",volmax);
   printf("\n");
   printf("======================================================================\n");
   printf("                  SOLVER AND FLOW RELATED OUTPUTS                     \n");
   printf("======================================================================\n");
   printf("\n");
   printf(" Mach number                %0.8e\n",s->mach);
   printf(" Reynolds number            %0.8e\n",s->rey);
   printf(" Alpha                      %0.8e\n",s->alpha);
   printf(" Beta                       %0.8e\n",s->beta);
   printf(" U infinity                 %0.8e\n",s->uinf);
   printf(" V infinity                 %0.8e\n",s->vinf);
   printf(" W infinity                 %0.8e\n",s->winf);
   printf(" Beta                       %0.8e\n",s->beta);
   printf(" CFL number                 %0.8e\n",s->cflnum);
   printf(" Number of steps            %d\n",s->nsteps);
   printf(" Time marching scheme       %s\n",s->scheme);
   // printf(" Time integration scheme    %s\n",g->timeInteg);
   printf(" Reconstruction accuracy    %d\n",g->order);

   printf("\n");
   printf("======================================================================\n");
}