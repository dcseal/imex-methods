This is a bunch of Matlab code that does stuff.

This is a crappy README file.



The driver for this whole shtuff is located in

    main.m           % main driver

User parameters and problem specific parameters are set here.  In order to
accomodate these paraemters, there's a global struct called params that gets
passed pretty much everywhere.



RK coefficients can be added by copying and modifying files coeffs* and
commenting the proper line in rk_integrator.m

User defined functions include:

    get_boundary.m   % Sets the ghost cells given a type of boundary condition
    
    qinit.m          % initial condition.
    auxinit.m        % initial condition for aux arrays if desired

    ConstructL.m     % This function depends on the function:
    fluxfunc.m       % function defining the problem




One of the most important parts of these IMEX routines is the matrix inversion
part.  If you create a new application, you need to modify the file

    rk_integrator.m  % main driver for taking a single time step

in order to accomodate the implicit solves.  Each implicit solve is of the
form:

    y = y_i + dt A_ii fI( tn + c[i], y );

Where the unknown variable is y.  The y_i come from the off diagonal part of
the butcher tableau, and fI is the implicit part of the system.


Whatever routine you use, it should should be able to handle any system of
ODES, provided the user specifies size and length of the array to be operated
on.  This routine depends on three functions: 

    coeffs_*         % butcher tableaus for additive rk scheme
    fE.m             % explicit part of the additive rk scheme
    fI.m             % implicit part of the additive rk scheme

The example provided in burger is an inviscous example where the implicit part
is actually linear.  In this case, one can get away with simply saving the
implicit matrix.  In more general cases, this needs to be modified.
