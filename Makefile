PROG =	pdnet

SRCS = fdriver.f90 pdnet.f90 pdnet_default.f90 pdnet_feasible.f90 pdnet_maxflow.f90\
	pdnet_read.f90 pdnet_solreport.f90 dblas1.f dnsubs.f

OBJS = fdriver.o pdnet.o pdnet_default.o pdnet_feasible.o pdnet_maxflow.o\
	pdnet_read.o pdnet_solreport.o dblas1.o dnsubs.o

LIBS =	

CC = gcc
CFLAGS = -C
#FC = ifort
FC = gfortran
#FFLAGS = -g -gline -u -dcfuns -nan -save -C -C=undefined
#FFLAGS = -check all -stand f95
FFLAGS = -O3
#F90 = ifort
F90 = gfortran
#F90FLAGS = -g -gline -u -dcfuns -nan -save -C -C=undefined
#F90FLAGS = -check all -stand f95
F90FLAGS = -O3
#LDFLAGS =  -g -gline -u -dcfuns -nan -save -C -C=undefined
#LDFLAGS = -check all -stand f95
LDFLAGS = -O3

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

