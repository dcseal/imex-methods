% Implicit right hand side function
function f = fI(t,q)


    global params
    mx   = params.mx;
    aux  = params.aux;
    meqn = params.meqn;

    q = reshape(q, mx, meqn );

    f = zeros( size(q) );
    f = reshape( -q/params.tau, mx*meqn, 1 );
    f = reshape(f, mx*meqn, 1 );


end
