# File definitions
SOURCES=sweeperUtils.f90 sweeper.f90 IO.f90 openmp.f90 openacc.f90 fspSolver.f90 EXACTISH.f90
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
	FLAGS=-std=f2003 -fall-intrinsics -Ofast -fbacktrace -fbounds-check -Wconversion -Wline-truncation 
else
	COMPILER=$(compiler) -c -o
	LINKER=$(compiler) -o
	FLAGS=-g -O0 -Mbounds
endif
# Object Definitions
OBJNAMES=$(SOURCES:.f90=.o)
OBJECTS=$(addprefix $(OBJDIR)/, $(OBJNAMES))
# Include Definitions
INCLUDE=-I$(OBJDIR)
# Executable Definition
EXE=EXACTISH.exe
# Command Definitions
RM=rm -rf
MK=mkdir

# Default target
.PHONY: 
all: $(OBJNAMES)
	$(LINKER) $(EXE) $(OBJECTS) $(FLAGS) $(INCLUDE)

# Object file target
%.o: %.f90 $(OBJDIR)
	$(COMPILER) $(OBJDIR)/$(@) $(FLAGS) $(INCLUDE) $< 

# Object directory target
$(OBJDIR):
	$(MK) $(OBJDIR)

# Clean target
.PHONY:
clean:
	$(RM) $(EXE)
	$(RM) $(OBJECTS)
	$(RM) $(OBJDIR)
