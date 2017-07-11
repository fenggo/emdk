# Makefile for GNU Linux / Debian stubs
#
SHELL = /bin/sh
#
# System-specific settings
#
#FC=ifort -static
#FC=gfortran -static -g -Wall -Wunused-parameter -fbounds-check 
FC=gfortran -static -g -Wall -ffpe-summary=underflow -fbounds-check  #none,all,underflow
CFLAG= -c 
#CFLAG= -c -fbacktrace -ftree-vectorize
#
# Link target
FOR=  mods.f90 \
      stretch.f90\
      analysis.f90 \
      memory.f90 \
      main.f90 supercell.f90 \
      control.f90 \
      write.f90 write_top.f90 write_uspex.f90 \
      write_moltop.f90 whether.f90 mmfit.f90\
      print_out.f90 plot.f90 \
      atom_type.f90 atom_type_hmx.f90 atom_type_cl20.f90 atom_type_nm.f90\
      project.f90 \
      neighbor.f90 \
      readin.f90 \
      string.f90 zerodata.f90 swing.f90\
      findmole.f90 putmol.f90 \
      gofr.f90 \
      zoom.f90 rotate.f90\
      defect.f90 \
      grid.f90 gulp_fit.f90 gulp_ff.f90  gulp_nw.f90\
      develop.f90 \
      move.f90 tinker.f90\
      velocity.f90 \
      dipole.f90 torsions.f90 torsion_angle.f90\
      dl_field.f90 angles.f90\
      structureDatabase.f90

OBJ= $(FOR:.f90=.o)
#
# Link rule
#
#
emdk: $(OBJ)
#	$(FC) $(OBJ)  -o exe/sym
#	$(FC) $(OBJ)  -static -o sym
	$(FC) $(OBJ)  -static -o emdk
#
#
# Compilation rules
%.o: %.f90
	$(FC) $(CFLAG) $<
#
#
.PHONY : clean
clean:
	rm -f *.o *.mod

