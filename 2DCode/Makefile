# Files
EXEC := Simulation
SRC  := $(wildcard *.c)
OBJ  := $(patsubst %.c, %.o, $(SRC))

# Options
CC 	:= gcc 
CFLAGS	:= -g -Wall -std=c99 
INCLUDE := -I/usr/include -I$$GSL_INC
LDFLAGS := -L/usr/lib -L$$GSL_LIB
LDLIBS 	:= -lgsl -lgslcblas  


# Rules
$(EXEC) : $(OBJ)
	$(CC) $(LDFLAGS) $(LDLIBS) -o $@ $^
%.o : %.c
	$(CC) $(CFLAGS) $(DEBUG) $(WD) $(INCLUDE) -c $<

main.o: SimulationSteps.h
extract_mesh_attributes.o: MeshAttributes.h mathfunctions.h
intialize.o: MeshAttributes.h
time_evolution.o: mathfunctions.h oneTimeStep.h
compute2DL.o: mathfunctions.h MeshAttributes.h
SlopeLimiter2D.o: mathfunctions.h MeshAttributes.h
boundary_conditions.o: MeshAttributes.h
extractXAxisEdges.o: MeshAttributes.h

# clean directive 'make clean'
.PHONY: clean neat echo
clean: neat
	$(RM) $(OBJ) $(EXEC)
neat:
	$(RM) *~.*~
