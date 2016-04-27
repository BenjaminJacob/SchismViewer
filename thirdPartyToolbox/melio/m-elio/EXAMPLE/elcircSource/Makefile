#
#  elcirc makefile
#  Mike Zulauf 04/30/03
#

#
# the executable and source code names (without extensions)
#
EXEC = ztest02g_da01c
ELCIRC = test02g_da01c
DSRC2 = dsrc2c
SFLUX = sflux_subs6

#
# comment out the following line if you don't want the build dependent
# upon the makefile, otherwise use the name of the makefile
#
# MAKEFILE = makefile

#
# location of HDF libraries
#
LIBDIR = /home/mazulauf/amb10xx/lib/x86

#
# compiler name
#
FC = ifc

#
# compilation flags
#
FFLAGS = -O3 -Bstatic

#
# the libraries
#
LIBS = -L$(LIBDIR) -lmfhdf -ldf -ljpeg -lz

#
# the object files
#
OBJS =  $(ELCIRC).o $(DSRC2).o $(SFLUX).o

#
# the actual build commands
#

$(EXEC): $(OBJS) $(MAKEFILE)
	$(FC) $(FFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

$(ELCIRC).o: $(ELCIRC).f90 $(MAKEFILE)
	$(FC) -c $(FFLAGS) $(ELCIRC).f90

$(DSRC2).o: $(DSRC2).f90 $(MAKEFILE)
	$(FC) -c $(FFLAGS) $(DSRC2).f90

$(SFLUX).o: $(SFLUX).f90 $(MAKEFILE)
	$(FC) -c $(FFLAGS) $(SFLUX).f90

#
# how to clean _all_ up
#

clean:
	rm -f *.o *.mod $(EXEC)
