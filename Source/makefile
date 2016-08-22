################################################################################
#Directory
WORK         =.
#
FC           = ifort
COMPILE.F    = $(FC) -c -traceback -check noarg_temp_created
LINK.F       = $(FC) -traceback
FFLAGS        = -O2 -check all
#
LIBS= -mkl -L$(MKLROOT)/lib 
#
#Objective functions
MYMAIN     = $(WORK)/mbi_main
#
#Currently used modules
MYGLOB     = mbi_glob
MYSUBFUN   = mbi_subfun
MYMISC     = mbi_misc
MYPRIOR    = mbi_prior
MYPOST     = mbi_post
#
#Bayesian related modules
GLOB       = $(WORK)/$(MYGLOB).o
MISC       = $(WORK)/$(MYMISC).o
SUBFUN     = $(WORK)/$(MYSUBFUN).o
PRIOR      = $(WORK)/$(MYPRIOR).o
POST       = $(WORK)/$(MYPOST).o
#
#Basic functions or subroutines
#
################################################################################
INCLUDE     = 
#
################################################################################
mainOBJ   = $(GLOB)  $(MISC)  $(SUBFUN)  $(PRIOR)  $(POST)  $(MYMAIN).o
#
################################################################################
clean: 
	rm -f $(mainOBJ) *.d *.o *.mod
#
################################################################################
mymain: $(INCLUDE) $(mainOBJ) 
	$(LINK.F) $(FFLAGS) $(mainOBJ) -o mymain $(LIBS) 
#
################################################################################
$(MYMAIN).o: $(INCLUDE) $(MYMAIN).f90
	$(COMPILE.F) $(FFLAGS) $(MYMAIN).f90 -o $@
$(PRIOR): $(INCLUDE) $(WORK)/$(MYPRIOR).f90
	$(COMPILE.F) $(FFLAGS) $(WORK)/$(MYPRIOR).f90 -o $@
$(POST): $(INCLUDE) $(WORK)/$(MYPOST).f90
	$(COMPILE.F) $(FFLAGS) $(WORK)/$(MYPOST).f90 -o $@
$(MISC): $(INCLUDE) $(WORK)/$(MYMISC).f90
	$(COMPILE.F) $(FFLAGS) $(WORK)/$(MYMISC).f90 -o $@
$(GLOB): $(INCLUDE) $(WORK)/$(MYGLOB).f90
	$(COMPILE.F) $(FFLAGS) $(WORK)/$(MYGLOB).f90 -o $@
$(SUBFUN): $(INCLUDE) $(WORK)/$(MYSUBFUN).f90
	$(COMPILE.F) $(FFLAGS) $(WORK)/$(MYSUBFUN).f90 -o $@
