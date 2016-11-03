# Files
EXEC := ChanNet 
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

main.o: ChannelsAndJunctions.h SimulationSteps.h
create_channels.o: ChannelsAndJunctions.h MeshAttributes.h
create_junctions.o: ChannelsAndJunctions.h MeshAttributes.h mathfunctions.h
intialize_channels.o: ChannelsAndJunctions.h mathfunctions.h MeshAttributes.h
initialize_junctions.o: ChannelsAndJunctions.h mathfunctions.h MeshAttributes.h
time_evolution.o: ChannelsAndJunctions.h mathfunctions.h oneTimeStep.h
computeL.o: ChannelsAndJunctions.h mathfunctions.h MeshAttributes.h
compute2DL.o: ChannelsAndJunctions.h mathfunctions.h MeshAttributes.h
minmod.o: ChannelsAndJunctions.h mathfunctions.h
SlopeLimiter2D.o: ChannelsAndJunctions.h mathfunctions.h
boundary_conditions.o: ChannelsAndJunctions.h MeshAttributes.h
constitutive_equations.o: mathfunctions.h

# clean directive 'make clean'
.PHONY: clean neat echo
clean: neat
	$(RM) $(OBJ) $(EXEC)
neat:
	$(RM) *~.*~
