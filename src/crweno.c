// ##############################################################
//
// crweno.c
//
// Compact WENO5 reconstruction scheme
// Written by Yong Su Jung
// ###############################################################
#include <stdio.h>
#include <stdlib.h>
#include "ham3dtypes.h"
#include "ham3dFunctionDefs.h"

void crweno(double **f,  //fluxes at each edge
            double **ql, // left states
            double **qr, // right states
            int    is,   // start index of chain
            int    ie,   // end index of chain
            double epsw,  // epsilon in limiter
            int    imax, // chainSize
            int    nq)   // number of variables
{
 
 int    i,n,im,k,j,imap;
 double qm2,qm1,qc,qp1,qp2,qp3;
 double f1,f2,f3,b1,b2,b3;
 double w1_opt,w2_opt,w3_opt;
 double a1,a2,a3,w1,w2,w3;
 double amat[imax-4],bmat[imax-4],cmat[imax-4],rmat[imax-4];
 double f0[imax];
  
  // for mapped weno (Henick et al)
  imap = 0;

  epsw = 1e-12;
  im = imax-2;
  
  //==========================
  // routine for left Biased 
  //===========================
  for (n=0;n<nq;n++) // number of flow variable
  {
     for (i=0;i<=ie;i++) // loop through the chain
     {   
       f0[i] = f[i][n];
     }
    
     //initialize to zero
     for(i=0;i<imax-4;i++)
     {
       amat[i] = 0.0;
       bmat[i] = 0.0;
       cmat[i] = 0.0;
       rmat[i] = 0.0;
     }

     for (i=is;i<im;i++)
     {

       qm2 = f0[i-2];
       qm1 = f0[i-1];
       qc  = f0[i];
       qp1 = f0[i+1];
       qp2 = f0[i+2];

       f1 = (qm1 + 5.0*qc)/6.0;
       f2 = (5.0*qc + qp1)/6.0;
       f3 = (qc + 5.0*qp1)/6.0;


     // smoothness indicators
       b1 = (13.0/12.0)*(qm2-2.0*qm1+qc)*(qm2-2.0*qm1+qc)+(1.0/4.0)*(qm2-4.0*qm1+3.0*qc)*(qm2-4.0*qm1+3.0*qc);
       b2 = (13.0/12.0)*(qm1-2.0*qc+qp1)*(qm1-2.0*qc+qp1)+(1.0/4.0)*(qm1-qp1)*(qm1-qp1);
       b3 = (13.0/12.0)*(qc-2.0*qp1+qp2)*(qc-2.0*qp1+qp2)+(1.0/4.0)*(3.0*qc-4.0*qp1+qp2)*(3.0*qc-4.0*qp1+qp2);


       // boudary face using WENO5
       if(i==is || i==im-1)
       {
         w1_opt = 0.1;
         w2_opt = 6.0/10.0;
         w3_opt = 3.0/10.0;

         a1 = w1_opt/pow((epsw+b1),2.0);
         a2 = w2_opt/pow((epsw+b2),2.0);
         a3 = w3_opt/pow((epsw+b3),2.0);

         w1 = a1/(a1+a2+a3);
         w2 = a2/(a1+a2+a3);
         w3 = a3/(a1+a2+a3);


         rmat[i-2] = w1/3.0*qm2-1.0/6.0*(7.0*w1+w2)*qm1+1.0/6.0*(11.0*w1+5.0*w2+2.0*w3)*qc + 1.0/6.0*(2.0*w2+5.0*w3)*qp1-w3/6.0*qp2; 
         amat[i-2] = 0.0;
         bmat[i-2] = 1.0;
         cmat[i-2] = 0.0;
       }
       else  // interior face using CRWENO5
       {
          w1_opt = 0.2;
          w2_opt = 0.5;
          w3_opt = 0.3;

          a1 = w1_opt/pow((epsw+b1),2.0);
          a2 = w2_opt/pow((epsw+b2),2.0);
          a3 = w3_opt/pow((epsw+b3),2.0);

          w1 = a1/(a1+a2+a3);
          w2 = a2/(a1+a2+a3);
          w3 = a3/(a1+a2+a3);

         if(imap==1) // mapping the weno weight
         {
           a1 = w1*(w1_opt+w1_opt*w1_opt -3.0*w1_opt*w1+w1*w1)/(w1_opt*w1_opt+w1*(1.0-2.0*w1_opt));
           a2 = w2*(w2_opt+w2_opt*w2_opt -3.0*w2_opt*w2+w2*w2)/(w2_opt*w2_opt+w2*(1.0-2.0*w2_opt));
           a3 = w3*(w3_opt+w3_opt*w3_opt -3.0*w3_opt*w3+w3*w3)/(w3_opt*w3_opt+w3*(1.0-2.0*w3_opt));
           // renomalize
           w1 = a1/(a1+a2+a3);
           w2 = a2/(a1+a2+a3);
           w3 = a3/(a1+a2+a3);
         }

          rmat[i-2] = w1*f1 + w2*f2 + w3*f3;
          amat[i-2] = (2.0/3.0)*w1 + (1.0/3.0)*w2;
          bmat[i-2] = (1.0/3.0)*w1 + (2.0/3.0)*(w2+w3);
          cmat[i-2] = (1.0/3.0)*w3;   
        }
      }

      //tridiagonal matrix solver
      tri(amat,bmat,cmat,rmat,imax-4);

      for (i=is;i<im;i++)
      { 
        ql[i][n] = rmat[i-2];
      }

   }

   //==========================
   // routine for right Biased 
   //===========================

   for (n=0;n<nq;n++) // number of flow variable
   {
     for (i=0;i<=ie;i++) // loop through the chain
     {   
       f0[i] = f[i][n];
     }

     for(i=0;i<imax-4;i++)
     {
       amat[i] = 0.0;
       bmat[i] = 0.0;
       cmat[i] = 0.0;
       rmat[i] = 0.0;
     }

     for (i=is;i<im;i++)
     {
       qm1 = f0[i-2];
       qc  = f0[i-1];
       qp1 = f0[i];
       qp2 = f0[i+1];
       qp3 = f0[i+2];

       f1 = (qp2 + 5.0*qp1)/6.0;
       f2 = (5.0*qp1 + qc)/6.0;
       f3 = (qp1 + 5.0*qc)/6.0;

        //smoothness indicators
        b1 = (13.0/12.0)*(qp3-2.0*qp2+qp1)*(qp3-2.0*qp2+qp1)+(1.0/4.0)*(qp3-4.0*qp2+3.0*qp1)*(qp3-4.0*qp2+3.0*qp1);
        b2 = (13.0/12.0)*(qp2-2.0*qp1+qc)*(qp2-2.0*qp1+qc)+(1.0/4.0)*(qp2-qc)*(qp2-qc);
        b3 = (13.0/12.0)*(qp1-2.0*qc+qm1)*(qp1-2.0*qc+qm1)+(1.0/4.0)*(3.0*qp1-4.0*qc+qm1)*(3.0*qp1-4.0*qc+qm1);

       // boudary face using WENO5
       if(i==is || i==im-1)
       {
         w1_opt = 0.1;
         w2_opt = 6.0/10.0;
         w3_opt = 3.0/10.0;

         a1 = w1_opt/pow((epsw+b1),2.0);
         a2 = w2_opt/pow((epsw+b2),2.0);
         a3 = w3_opt/pow((epsw+b3),2.0);

         w1 = a1/(a1+a2+a3);
         w2 = a2/(a1+a2+a3);
         w3 = a3/(a1+a2+a3);

         rmat[i-2] = w1/3.0*qp3-1.0/6.0*(7.0*w1+w2)*qp2+1.0/6.0*(11.0*w1+5.0*w2+2.0*w3)*qp1 + 1.0/6.0*(2.0*w2+5.0*w3)*qc-w3/6.0*qm1; 
         amat[i-2] = 0.0;
         bmat[i-2] = 1.0;
         cmat[i-2] = 0.0;
       }
       else  // interior face
       {
         //optinal weight
         w1_opt = 0.2;
         w2_opt = 0.5;
         w3_opt = 0.3;
         
         //limiting
         a1 = w1_opt/pow((epsw+b1),2.0);
         a2 = w2_opt/pow((epsw+b2),2.0);
         a3 = w3_opt/pow((epsw+b3),2.0);

         //making the weights convex
         w1 = a1 / (a1+a2+a3);
         w2 = a2 / (a1+a2+a3);
         w3 = a3 / (a1+a2+a3);

         if(imap==1) // mapping the weno weight
         {
           a1 = w1*(w1_opt+w1_opt*w1_opt -3.0*w1_opt*w1+w1*w1)/(w1_opt*w1_opt+w1*(1.0-2.0*w1_opt));
           a2 = w2*(w2_opt+w2_opt*w2_opt -3.0*w2_opt*w2+w2*w2)/(w2_opt*w2_opt+w2*(1.0-2.0*w2_opt));
           a3 = w3*(w3_opt+w3_opt*w3_opt -3.0*w3_opt*w3+w3*w3)/(w3_opt*w3_opt+w3*(1.0-2.0*w3_opt));
           // renomalize
           w1 = a1/(a1+a2+a3);
           w2 = a2/(a1+a2+a3);
           w3 = a3/(a1+a2+a3);
         }


         rmat[i-2] = w1*f1 + w2*f2 + w3*f3;
         cmat[i-2] = (2.0/3.0)*w1 + (1.0/3.0)*w2;
         bmat[i-2] = (1.0/3.0)*w1 + (2.0/3.0)*(w2+w3);
         amat[i-2] = (1.0/3.0)*w3;
       }

     }
     
     //tridiagonal matrix solver
     tri(amat,bmat,cmat,rmat,imax-4);

     // n is for the  variables
     for (i=is;i<im;i++)
     { 
       qr[i][n] = rmat[i-2];
     }
     
 } //loop for number of flow variable

} //end function


