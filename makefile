# File definitions
SOURCES=sweeperUtils.f90 sweeper.f90 IO.f90 openmp.f90 openacc.f90 fspSolver.f90 
# Directory Definitions
EMPTY=
ifeq ($(source),$(EMPTY))
	VPATH=/nfs-home/aarograh/homework/eecs587/EXACT-ISH
else
	VPATH=$(source)
endif
ifeq ($(build),$(EMPTY))
	OBJDIR=$(PWD)/obj
else
	OBJDIR=$(build)
endif
# Compiler Stuff
ifeq ($(compiler),$(EMPTY))
	COMPILER=gfortran -J$(OBJDIR) -c -o
	LINKER=gfortran -o
#	COMPILER=/nfs-home/aarograh/homework/eecs587/linux86-64/14.10/bin/pgfortran -c -o
#	LINKER=/nfs-home/aarograh/homework/eecs587/linux86-64/14.10/bin/pgfortran -o
else
	COMPILER=$(compiler)
endif
FLAGS=-std=f2003 -Wall -fall-intrinsics -Ofast -g -fbacktrace -fbounds-check -Wall
# Object Definitions
OBJNAMES=$(SOURCES:.f90=.o)
OBJECTS=$(addprefix $(OBJDIR)/, $(OBJNAMES))
JIPUOBJ=jipu_main.o
AARONOBJ=aaron_main.o
# Include Definitions
INCLUDE=-I$(OBJDIR)
# Executable Definition
EXECUTABLES=aaron jipu
EXE1=jipu_main.exe $(OBJDIR)/jipu_main.o
EXE2=aaron_main.exe $(OBJDIR)/aaron_main.o
# Command Definitions
RM=rm -rf
MK=mkdir

# Default target
.PHONY: 
all: $(EXECUTABLES)

jipu: $(OBJNAMES) $(JIPUOBJ)
	$(LINKER) $(EXE1) $(FLAGS) $(INCLUDE) $(OBJECTS)

aaron: $(OBJNAMES) $(AARONOBJ)
	$(LINKER) $(EXE2) $(FLAGS) $(INCLUDE) $(OBJECTS)

# Object file target
%.o: %.f90 $(OBJDIR)
	$(COMPILER) $(OBJDIR)/$(@) $(FLAGS) $(INCLUDE) $< 

# Object directory target
$(OBJDIR):
	$(MK) $(OBJDIR)

# Clean target
.PHONY:
clean:
	$(RM) $(EXE2)
	$(RM) $(OBJECTS)
	$(RM) $(OBJDIR)
