# File definitions
SOURCENAMES=sweeper.f90 
MAIN=$(SOURCEDIR)/EXACT-ISH.f90
# Compiler Stuff
COMPILER=gfortran
FLAGS=-std=f2003 -Wall -fall-intrinsics -Ofast -g -fbacktrace -fbounds-check -Wall
# Directory Definitions
BUILDDIR=${PWD}
SOURCEDIR=/nfs-home/aarograh/homework/eecs587/EXACT-ISH
OBJDIR=$(BUILDDIR)/obj
# Source Definitions
SOURCES=$(addprefix $(SOURCEDIR)/, $(SOURCENAMES))
# Object Definitions
OBJNAMES=$(SOURCENAMES:.f90=.o)
OBJECTS=$(addprefix $(OBJDIR)/, $(OBJNAMES))
# Include Definitions
INCLUDE=-I$(OBJDIR)
# Executable Definition
EXECUTABLE=EXACT-ISH.exe
# Command Definitions
MV=mv
RM=rm -rf
MK=mkdir

# Default target
all: $(OBJDIR) $(OBJECTS)
	$(COMPILER) $(FLAGS) $(OBJECTS) $(MAIN) $(INCLUDE) -o $(EXECUTABLE)

# Object file target
$(OBJECTS): %.o : $(SOURCES)
	$(COMPILER) $(FLAGS) $< -c; $(MV) *.o $(OBJDIR); $(MV) *.mod $(OBJDIR)

# Object directory target
$(OBJDIR):
	$(MK) $(OBJDIR)

# Clean target
clean:
	$(RM) EXACT-ISH.exe; $(RM) $(OBJDIR)
