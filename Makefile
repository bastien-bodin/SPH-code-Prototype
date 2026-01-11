# Compiler settings
FC = gfortran

# Base flags (Common to both versions)
# -O3: High optimization
# -march=native: Optimize for the current CPU
FFLAGS_BASE = -O3 -march=native -Wall

# --- OpenMP Toggle Logic ---
# Usage: 
#   make OMP=1 (default, parallel)
#   make OMP=0 (serial)
OMP ?= 1

ifeq ($(OMP), 1)
    FFLAGS = $(FFLAGS_BASE) -fopenmp
    msg = "Building PARALLEL version (OpenMP enabled)"
else
    FFLAGS = $(FFLAGS_BASE)
    msg = "Building SERIAL version (OpenMP disabled)"
endif

# Source files
SRCS = parameters.f90 \
       particles.f90 \
       kernels.f90 \
       sort_parts.f90 \
       get_neighbours.f90 \
       equations.f90 \
       density.f90 \
       forces.f90 \
       integrator.f90 \
       geometries.f90 \
       application.f90 \
       main.f90

OBJS = $(SRCS:.f90=.o)
TARGET = sph_europa

all: info $(TARGET)

info:
	@echo "------------------------------------------"
	@echo $(msg)
	@echo "FFLAGS: $(FFLAGS)"
	@echo "------------------------------------------"

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# Dependencies (Order matters for .mod files)
particles.o: parameters.o
kernels.o: parameters.o
sort_parts.o: parameters.o particles.o
get_neighbours.o: parameters.o particles.o sort_parts.o
equations.o: parameters.o particles.o
density.o: parameters.o particles.o kernels.o
forces.o: parameters.o particles.o kernels.o equations.o
integrator.o: parameters.o particles.o
geometries.o: parameters.o particles.o
application.o: parameters.o particles.o sort_parts.o get_neighbours.o \
               kernels.o equations.o density.o forces.o integrator.o geometries.o
main.o: application.o

clean:
	rm -f $(OBJS) *.mod $(TARGET) output_*.csv