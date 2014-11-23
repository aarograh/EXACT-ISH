# File definitions
SOURCES=sweeperUtils.f90 sweeper.f90 fspSolver.f90 IO.f90 EXACT-ISH.f90
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
COMPILER=gfortran -J$(OBJDIR) -c -o
LINKER=gfortran -o
FLAGS=-std=f2003 -Wall -fall-intrinsics -Ofast -g -fbacktrace -fbounds-check -Wall
# Object Definitions
OBJNAMES=$(SOURCES:.f90=.o)
OBJECTS=$(addprefix $(OBJDIR)/, $(OBJNAMES))
# Include Definitions
INCLUDE=-I$(OBJDIR)
# Executable Definition
EXECUTABLE=EXACT-ISH.exe
# Command Definitions
RM=rm -rf
MK=mkdir

# Default target
.PHONY: 
all: $(OBJNAMES) 
	@$(LINKER) $(EXECUTABLE) $(FLAGS) $(INCLUDE) $(OBJECTS)

# Object file target
%.o: %.f90 $(OBJDIR)
	@$(COMPILER) $(OBJDIR)/$(@) $(FLAGS) $(INCLUDE) $< 

# Object directory target
$(OBJDIR):
	@$(MK) $(OBJDIR)

# Clean target
.PHONY:
clean:
	@$(RM) $(EXECUTABLE)
	@$(RM) $(OBJECTS)
	@$(RM) $(OBJDIR)
