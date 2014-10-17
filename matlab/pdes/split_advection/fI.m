% Implicit right hand side function
function f = fI(t,q)


    global params
    A    = params.A;
    mx   = params.mx;
    meqn = params.meqn;

    q = reshape(q, mx, meqn );

    f = zeros( size(q) );
    f = params.uimp * ConstructL( q, t );
    f = reshape(f, mx*meqn, 1 );

end
