include ../Makefile.config

CXXFLAGS_EXTRA	?= -I.

SRCS              = \
    integrate.cpp \
    tensors.cpp   \
    set_boundary.cpp \
    construct_l.cpp \
    global_main.cpp \
    req_state.cpp \
    rk_integrator.cpp \
    create_output_dir.cpp \
    solve.cpp \
    output.cpp

OBJS              = $(SRCS:.cpp=.o)

lib:              all

default:          all

all:              lib$(LIBNAME).a

lib$(LIBNAME).a:  $(OBJS) 
		  ${AR} -r $@ $(OBJS)

#### This part doesn't work ... ###
#all:              lib$(LIBNAME).so

#lib$(LIBNAME).so:  $(OBJS) 
#		  ${LD} -G -fPIC $@ $(OBJS) -o lib$(LIBNAME).so
####################################


clean:          
		  $(RM) *.o *.a *.so

