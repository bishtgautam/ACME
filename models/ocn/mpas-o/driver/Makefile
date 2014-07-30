.SUFFIXES: .F .o

OBJS = ocn_comp_mct.o \
       ocn_communicator.o \
       mpaso_cpl_indices.o \
	   mpaso_mct_vars.o

all: $(OBJS)

ocn_comp_mct.o: ocn_communicator.o mpaso_cpl_indices.o mpaso_mct_vars.o

ocn_communicator.o:

mpaso_cpl_indices.o:

mpaso_mct_vars.o:

clean:
	$(RM) *.o *.mod *.f90

.F.o:
	$(FC) -c $(CPPFLAGS) $(CPPINCLUDES) $(FFLAGS) $(FCINCLUDES) $(FREEFLAGS)  -I../framework -I../core_$(CORE) -I../external/esmf_time_f90 $<

#.F.o:
#	$(RM) $@ $*.mod
#	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
#	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../core_$(CORE) -I../external/esmf_time_f90
