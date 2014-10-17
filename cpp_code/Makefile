include ./Makefile.config

# add all the directories here.  lib is the one everyone cares about!
DIRS = \
    lib \
    advection \
    burgers \
    stiff_advection_reaction \
    test_integrate \
    diffusion

default:	what

all:		
		for i in $(DIRS); do (cd $$i && $(MAKE) $@) || exit 1; done

# this allows you to do something like "make lib"
.PHONY:		$(DIRS)
$(DIRS):
		(cd $@ && $(MAKE)) || exit 1

clean:
		for i in $(DIRS); do (cd $$i && $(MAKE) $@) || exit 1; done

what:
		@echo "all lib clean"
