#############################################################
#                                                           #
#    Original code at Github: xiaodongli1986/LSS_Code       #
#         hacked to calculate just Beta-Skeleton            #
#                                                           #
#############################################################


EXE2 =  ../bin/LSS_BSK_calc

OBJS = LSS_tools.o  LSS_cosmo_funs.o LSS_settings_init.o LSS_cell.o LSS_smooth.o LSS_BSK.o
F90C = gfortran
F90FLAGS = 

default: $(EXE2)
all: $(EXE2)
bsk: $(EXE2)
BSK: $(EXE2)


$(EXE2): $(OBJS) LSS_main_BSK_calc.o
	$(F90C) -o $(EXE2) $(OBJS) LSS_main_BSK_calc.o

LSS_cosmo_funs.o: LSS_tools.o
LSS_settings_init.o: LSS_cosmo_funs.o 
LSS_cell.o: LSS_settings_init.o
LSS_smooth.o: LSS_cell.o
LSS_grad_fields.o: LSS_smooth.o LSS_mu_tools.o
LSS_chisq.o: LSS_grad_fields.o 
LSS_BSK.o : LSS_chisq.o
LSS_main_BSK_calc.o: LSS_BSK.o 


%.o: %.f90
	$(F90C) -c $*.f90 $(F90FLAGS)
%.o: %.F90
	$(F90C) -c $< -o $@

clean :
	rm -rf *.o *.mod $(EXE2) 
