#####################################################################
#                                                                   #
#                       MAKEFILE FOR HAMSTR 3D                      #
#                                                                   #
#####################################################################
#F90 = ifort
F90 = gfortran
#CC  = icc
CC  = gcc
MPIFC = /usr/local/bin/mpicc

#FFLAGS =  -w -O2 -r8 
FFLAGS  =  -w -O2 -fdefault-real-8 # uncomment this for gfortran
CFLAGS =  -w -O2 #-g #-traceback  -check uninit  #-warn-all#-g #-O2 #-g

#FFLAGS =  -w -O2 -lgfortran -fdefault-real-8 # uncomment this for gfortran
#CFLAGS =  -w -O2 -fbacktrace -fbounds-check # -g #-traceback  -check uninit  #-w


LFLAGS =$(CFLAGS) -lm
SRCDIR = ./src
OBJECTS = $(SRCDIR)/ham3d.o \
          $(SRCDIR)/readGrid.o \
          $(SRCDIR)/preprocess.o \
          $(SRCDIR)/find_faces.o\
	  		 $(SRCDIR)/initflow.o \
	  		 $(SRCDIR)/stepSolution.o \
	  		 $(SRCDIR)/updateSoln.o \
	  		 $(SRCDIR)/communication.o\
	       $(SRCDIR)/communicationLinear.o \
	       $(SRCDIR)/computeRHS.o \
	       $(SRCDIR)/computeRHSv.o\
          $(SRCDIR)/computeRHSk.o \
          $(SRCDIR)/computeRHSkv.o \
          $(SRCDIR)/computeLinearRHS.o \
          $(SRCDIR)/flux_roe3d.o \
          $(SRCDIR)/wallFlux.o\
          $(SRCDIR)/flux_roe2d.o \
          $(SRCDIR)/periodic_bc.o \
          $(SRCDIR)/apply_periodic.o \
          $(SRCDIR)/apply_periodic_LHS.o\
	       $(SRCDIR)/muscld.o \
	       $(SRCDIR)/muscld_deriv.o \
	       $(SRCDIR)/roeflx.o \
	       $(SRCDIR)/computeForce.o \
	       $(SRCDIR)/outputSolution.o \
	       $(SRCDIR)/smoothGrid.o \
	       $(SRCDIR)/jac_roe.o \
	       $(SRCDIR)/flux_visc.o \
	       $(SRCDIR)/triSolvers.o \
	       $(SRCDIR)/mathops.o \
	       $(SRCDIR)/findDiagonals.o\
	       $(SRCDIR)/ADI.o \
	       $(SRCDIR)/DADI.o \
	       $(SRCDIR)/jac_visc.o \
	       $(SRCDIR)/gaussSeidel.o \
	       $(SRCDIR)/lineGaussSeidel.o \
	       $(SRCDIR)/weno.o \
	       $(SRCDIR)/weno_deriv.o \
	       $(SRCDIR)/crweno.o \
	       $(SRCDIR)/tri.o

# Link Instruction
ham2d: $(OBJECTS)
	$(MPIFC) $(OBJECTS) $(LFLAGS) -o ham3d

clean:
	@rm -rf *.o *.mod *.*~ ham3d

%.o:%.F90
	$(F90) $(FFLAGS) -c $< -o $*.o

%.o:%.f90
	$(F90) $(FFLAGS) -c $< -o $*.o
%.o:%.c
	$(MPIFC) $(CFLAGS) -c $< -o $*.o

%.o.:%.mod	

#####################################################################
#                         END OF FILE                               #
#####################################################################
