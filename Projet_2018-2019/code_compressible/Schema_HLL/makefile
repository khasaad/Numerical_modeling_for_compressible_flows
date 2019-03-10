FORTRAN=gfortran
#COPTS= -c  -O -pedantic -fbounds-check -Wall -g -fbacktrace 
COPTS= -c  -O -pedantic -fbounds-check 
#-c -checkbound -O -r8 -v
PROG=run
SUFFIXES=.f90.o
.SUFFIXES : .f90 .o
$(SUFFIXES):
	$(FORTRAN) $(COPTS) $*.f90

OBJS	      = fonction_HLL.o\
		main_HLL.o


$(PROG): $(OBJS)
	$(FORTRAN) $(OBJS) -o $(PROG)

clean:;	rm -f *.o *.mod *~ run resul* density.pdf energy.pdf velocity.pdf



