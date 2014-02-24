

default: all

mission:
	f2py -c --f90flags='$(FFLAGS)' -m mission mission.f90

all:
	f2py -c --f90flags='$(FFLAGS)' -m mission mission.f90

clean:
	-rm *.o
	-rm *.mod
