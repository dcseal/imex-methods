include ../Makefile.config

# gnumake thing: this appends it if it isn't defined (?)
# this either defines it or appends it if it isn't already defined
CXXFLAGS_EXTRA	?= -I../lib
LDFLAGS_EXTRA	?= -L../lib -L$(PETSC_DIR)/lib

#LDLIBS_EXTRA	?= -l$(LIBNAME)
LDLIBS_EXTRA	?= -l$(LIBNAME) -lpetscksp -lpthread -lrt -lmpich -ldl -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc

SRCS		= \
    main.cpp \
    init_params.cpp \
    init_state.cpp \
    user_funcs.cpp

OBJS		= $(SRCS:.cpp=.o)
EXE		= main

default:	all

all:		$(EXE)

$(EXE):		lib $(OBJS) 
		$(CXX) $(LDFLAGS) $(OBJS) $(LDLIBS) -o $(EXE)

lib:
		(cd ../$@ && $(MAKE) $@ && cd ..) || exit 1;

clean:          
		$(RM) *.o $(EXE)
