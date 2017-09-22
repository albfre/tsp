CC = g++
CFLAGS = -g -funroll-loops -std=gnu++14 -O2 -Wall -pedantic
SRCS = \
       TravelingSalespersonProblemSolver.cpp \
       TravelingSalespersonProblemSolver.h \
       Distances.cpp \
       Distances.h \

OBJS = TravelingSalespersonProblemSolver.o Distances.o
TEST = TestTSP
TESTSRC = TestTSP.cpp

TravelingSalespersonProblemSolver.cpp : InitialTours.h HelperFunctions.h NeighborAlgorithms.h PerformHelsgaunMove.h

all: product test

product: $(SRCS) \
; $(CC) $(CFLAGS) -c $(SRCS)

test: $(TESTSRC) \
; $(CC) $(CFLAGS) -o $(TEST) $(TESTSRC) $(OBJS)

clean: \
; $(RM) $(OBJS)

write: \
; echo $(OBJS)
