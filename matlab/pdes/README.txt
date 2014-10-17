The structure of this directory is very similar to that of matlab/odes
directory.

There is one library 'matlab/pdes/lib/' which contains all of the common library
routines.  Each application is in its own folder and calls these library
routines.

The file: main.m is the main driver for each application.  This section is not
quite as clean as the ode section beacuse the applications use slightly
different flags that haven't been cleaned up.

Parameters are again set in a file called set_params.m

A few 2D examples have been written.  They can be found in matlab/pdes/2d/,
and that directory also follows the same file structure.
