The matlab applications are broken up into an ODE section and a PDE section.
The library applications are located in odes/lib and pdes/lib, and respective
applications can be found in each of the above mentioned folders.

Each applicaiton in the PDES section is expected to provide the following
functions:

    fE.m, fI.m, implicit_solve.m, main.m, auxinit.m, qinit.m, qexact.m

as well as the parameters file:

    set_params.m.

The format for hte function implicit_solve.m is the following:


    function y = implicit_solve( ys, ts, as )
        % Function providing the implicit solve for:
        % y = ys + as * fI( y, t );

        % Do some stuff here.  If there is not implicit part to the problem, the
        % proper way to initialize this is by:
        y = yexp;

    end


