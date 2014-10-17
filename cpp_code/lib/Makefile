include ../Makefile.config

# CXXFLAGS_EXTRA	?= -I. -I$(PETSC_DIR)/include -D IMPLICIT
CXXFLAGS_EXTRA	?= -I. -I$(PETSC_DIR)/include

SRCS              = \
    integrate.cpp \
    tensors.cpp   \
    set_boundary.cpp \
    construct_l.cpp \
    global_main.cpp \
    req_state.cpp \
    rk_integrator.cpp \
    sdc_integrator.cpp \
    ResCoeffs.cpp \
    create_output_dir.cpp \
    solve.cpp \
    matrix_solve.cpp \
    output.cpp

OBJS              = $(SRCS:.cpp=.o)
LDFLAGS_EXTRA     = $(PETSC_DIR)/lib/

lib:              all

default:          all

all:              lib$(LIBNAME).a

lib$(LIBNAME).a:  $(OBJS) 
		  ${AR} -r $@ $(OBJS)

###################################
#### This part doesn't work ... ###
#all:              lib$(LIBNAME).so

#lib$(LIBNAME).so:  $(OBJS) 
#		  ${LD} -G -fPIC $@ $(OBJS) -o lib$(LIBNAME).so
####################################


clean:          
		  $(RM) *.o *.a *.so

