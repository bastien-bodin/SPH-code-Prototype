#compiler
FC = gfortran

# compile flags
FCFLAGS = -g -c -fdefault-real-8 -fbacktrace -fno-align-commons -O3
# link flags
FLFLAGS = -g -fbacktrace


# edition de liens
main: main.o parameters.o particles.o kernels.o sort_parts.o equations.o get_neighbours.o application.o
	$(FC) -o main.x *.o
# compilation
main.o: main.f90 application.o
	$(FC) $(FCFLAGS) main.f90
application.o: application.f90 get_neighbours.o sort_parts.o equations.o kernels.o particles.o parameters.o
	$(FC) $(FCFLAGS) application.f90
get_neighbours.o: get_neighbours.f90 sort_parts.o equations.o particles.o parameters.o
	$(FC) $(FCFLAGS) get_neighbours.f90
equations.o: equations.f90 particles.o parameters.o
	$(FC) $(FCFLAGS) equations.f90
sort_parts.o: sort_parts.f90 particles.o parameters.o
	$(FC) $(FCFLAGS) sort_parts.f90
kernels.o: kernels.f90 parameters.o
	$(FC) $(FCFLAGS) kernels.f90
particles.o: particles.f90 parameters.o
	$(FC) $(FCFLAGS) particles.f90
parameters.o: parameters.f90
	$(FC) $(FCFLAGS) parameters.f90

clean:
	rm *.o *.mod *.x
