% Explicit right hand side function
function f = fE(t,q)

    global params
    mx   = params.mx;
    meqn = params.meqn;

    % source term gets treated implicitly.
    q = reshape( q, mx, meqn );
    f = zeros( size(q) );
    f = params.uexp * ConstructL( q, t );
    f = reshape( f, mx*meqn, 1 );

end
