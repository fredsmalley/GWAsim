## Makefile -- GLN makefile
##
## Joe Song
## Created: September 27, 2008
## Modified: 
##   February 26, 2010
##   October 4, 2011
##   November 2, 2011
##   March 11, 2012.  
##      Added support for C++11 
##   March 13, 2012. 
##      - Changed compilation on MacOS using Apple LLVM Compiler 
##        instead of GCC
##      - Corrected options for C++11 support
##   May 1, 2013. Older version of the CGLN files are removed from the 
##      Makefile.  They can still be retrieved using the time stamp of 
##      April 30, 2013

PROGMPI = glnmpi
PROGSP = glnsp
PROGCGLN = cgln

UNAME := $(shell uname -s)

ifeq ($(MPI), mpi)
	# Linux cluster
	CC = mpicc          
	CXX = mpicxx
	MPIFLAG = -D ENABLE_MPI
	PROG = $(PROGMPI)
else ifeq ($(MPI), openmpi)
	# Mac OS Lion OpenMPI
	CC = openmpicc  
	CXX = openmpicxx
	PROG = $(PROGMPI)
	MPIFLAG = -D ENABLE_MPI
else ifndef $(MPI)
ifeq ($(UNAME), Darwin)
	CC = clang               
	CXX = clang++ # to support C++11
else
	CC = gcc
	CXX = g++
endif
	MPIFLAG =
	PROG = $(PROGSP)
endif

# C++11 Setup:
ifeq ($(UNAME), Darwin) # on MacOS
    C++11 := -std=c++0x -stdlib=libc++
    C++11LIB := -lc++
    BOOSTPATH := -I /usr/local/include
else # on other OS (Linux)
    C++11 := -std=c++0x
    C++11LIB := 
    BOOSTPATH := -I /usr/local/include
    BOOSTLIB := 
endif

.SUFFIXES:
.SUFFIXES: .c .cpp .o .h

RM = rm -f
STRIP = strip -v

#EXCEPTION = -fhandle-exceptions
EXCEPTION =

# Debug flags
DFLAG = -Wall $(MPIFLAG)
#DFLAG = -g -Wall
#DFLAG = -DDebug1 -DDebug2 -DDebug3

# Includes, libraries and ldflags etc.
INCS =
##INCS =  -I/usr/include

#CXXFLAGS = -O3 $(C++11) $(DFLAG) $(INCS) $(BOOSTPATH) $(BOOSTLIB) $(EXCEPTION)
CXXFLAGS = -O3 $(C++11) $(DFLAG) $(INCS) $(EXCEPTION) $(BOOSTPATH)


#LIB = -L /opt/local/lib
LIB = $(C++11LIB)

LDFLAGS =

# Object and Source files for Generalized Logical Network Modeling
SRCSCPP = \
    AdjustPvalue.cpp \
    binom.cpp \
    ChisqDirTest.cpp ChisqTests.cpp ComparativeChisq.cpp CompareChisq.cpp \
    ComparePathways.cpp ConditionalAssociation.cpp \
    effective_dim.cpp EMTFisherTests.cpp EnumMethod.cpp \
    EnvGLNCmp.cpp EnvGLNCmpIO.cpp EnvGLNCmpParentSets.cpp \
    EnvGLNRec.cpp EnvGLNRecIO.cpp EnvGLNRecParentSets.cpp \
    ExactCPChisq2.cpp ExactCPChisq3.cpp ExactCPChisqTesting.cpp \
    ExactMultiTableTest.cpp \
    fexact.cpp \
    GeneralizedTruthTable.cpp GLN.cpp GLN_evaluate.cpp GLN_generate.cpp \
    GLN_IO.cpp GLN_manpage.cpp GLN_simulate.cpp GLN_Test.cpp \
    GLNAlcohol.cpp GLNCompare.cpp GLNPureReconstruction.cpp \
    GLNReconstruct.cpp \
    InteractionEvaluationController.cpp \
    JobGLN.cpp JobResultCmp.cpp JobResultRec.cpp \
    main.cpp \
    Node.cpp \
    PathwayStats.cpp PermStatTable.cpp \
    ResampleGLN.cpp ResampleTCF.cpp RunMPIJobs.cpp RunSPJobs.cpp \
    StatDistributions.cpp SetOps.cpp sum_chisqs.cpp \
    Topology.cpp TrajectoryCollection.cpp \
    TrajectoryCollectionIO.cpp Transition.cpp TransitionTable.cpp \
    TransitionTableIO.cpp \
    ExactFunctionalTest.cpp

OBJSCPP = \
    AdjustPvalue.o \
    binom.o \
    ChisqDirTest.o ChisqTests.o ComparativeChisq.o CompareChisq.o \
    ComparePathways.o ConditionalAssociation.o \
    effective_dim.o EMTFisherTests.o EnumMethod.o \
    EnvGLNCmp.o EnvGLNCmpIO.o EnvGLNCmpParentSets.o \
    EnvGLNRec.o EnvGLNRecIO.o EnvGLNRecParentSets.o \
    ExactCPChisq2.o ExactCPChisq3.o ExactCPChisqTesting.o \
    ExactMultiTableTest.o \
    fexact.o \
    GeneralizedTruthTable.o GLN.o GLN_evaluate.o GLN_generate.o \
    GLN_IO.o GLN_manpage.o GLN_simulate.o GLN_Test.o \
    GLNAlcohol.o GLNCompare.o GLNPureReconstruction.o \
    GLNReconstruct.o \
    InteractionEvaluationController.o \
    JobGLN.o JobResultCmp.o JobResultRec.o \
    main.o \
    Node.o \
    PathwayStats.o PermStatTable.o \
    ResampleGLN.o ResampleTCF.o RunMPIJobs.o RunSPJobs.o \
    StatDistributions.o SetOps.o sum_chisqs.o \
    Topology.o TrajectoryCollection.o \
    TrajectoryCollectionIO.o Transition.o TransitionTable.o \
    TransitionTableIO.o \
    ExactFunctionalTest.o

CGLNCPPS = AdjustPvalue.cpp ChisqTests.cpp ChisqDirTest.cpp CompareChisq.cpp \
    DS_AdjacencyList.cpp\
    DS_CompulsiveParents.cpp \
    DS_LocalConstraint.cpp DS_ConstraintFile.cpp\
    DS_main.cpp  effective_dim.cpp \
    EffCGLN_GlobalConstraint.cpp EffCGLN_DynamicPattern.cpp \
    EMTFisherTests.cpp EnvGLNRec.cpp EnvGLNRecIO.cpp EnvGLNRecParentSets.cpp \
    EnumMethod.cpp ExactMultiTableTest.cpp fexact.cpp \
    GLNPureReconstruction.cpp \
    GeneralizedTruthTable.cpp GLN.cpp GLN_IO.cpp GLN_simulate.cpp \
    GLNReconstruct.cpp \
    JobGLN.cpp JobResultRec.cpp \
    kStop.cpp kStop_Node.cpp kStop_SegmentPath.cpp kStop_TestCase.cpp\
    Node.cpp PathwayStats.cpp \
    PermStatTable.cpp \
    RunSPJobs.cpp RunMPIJobs.o \
    StatDistributions.cpp SetOps.cpp sum_chisqs.cpp \
    Topology.cpp \
    Transition.cpp TransitionTable.cpp \
    TrajectoryCollection.cpp TrajectoryCollectionIO.cpp \
    InteractionEvaluationController.cpp ConditionalAssociation.cpp

CGLNOBJS = AdjustPvalue.o ChisqTests.o ChisqDirTest.o CompareChisq.o \
    DS_AdjacencyList.o\
    DS_CompulsiveParents.o \
    DS_LocalConstraint.o DS_ConstraintFile.o\
    DS_main.o  effective_dim.o \
    EffCGLN_GlobalConstraint.o EffCGLN_DynamicPattern.o \
    EMTFisherTests.o EnvGLNRec.o EnvGLNRecIO.o EnvGLNRecParentSets.o \
    EnumMethod.o ExactMultiTableTest.o fexact.o \
    GLNPureReconstruction.o \
    GeneralizedTruthTable.o GLN.o GLN_IO.o GLN_simulate.o \
    GLNReconstruct.o \
    JobGLN.o JobResultRec.o \
    kStop.o kStop_Node.o kStop_SegmentPath.o kStop_TestCase.o\
    Node.o PathwayStats.o \
    PermStatTable.o\
    RunSPJobs.o RunMPIJobs.o \
    StatDistributions.o  SetOps.o sum_chisqs.o \
    Topology.o \
    Transition.o TransitionTable.o \
    TrajectoryCollection.o TrajectoryCollectionIO.o \
    InteractionEvaluationController.o ConditionalAssociation.o


#
# make the various executables
#

PROGRAM = $(PROG)

all :: $(PROGRAM) $(PROGCGLN)

touch:
	$@ $(SRCSCPP) $(CGLNCPPS)

$(PROGCGLN): $(CGLNOBJS)
	- $(RM) $@
	$(CXX) -o $@ $(CGLNOBJS) $(LIB) $(LDFLAGS)

$(PROGRAM): $(OBJSCPP)
	- $(RM) $@
	$(CXX) -o $@ $(OBJSCPP) $(LIB) $(LDFLAGS)

%.o : %.cpp %.h
	$(CXX) -c $(CXXFLAGS) $< -o $@

##%.o : %.c %.h
##	$(CC) -c $(CXXFLAGS) $< -o $@

strip:
	$(STRIP) $(PROGRAM) $(PROGCGLN)

.PHONY : clean
clean:
	- $(RM) *.o *.log $(PROGSP) $(PROGMPI) $(PROGCGLN) core

depend:
	makedepend -- $(CXXFLAGS) $(INCS) -- $(SRCSCPP) $(SRCSC) $(CGLNCPPS)

