% Explicit right hand side function
function f = fE(t,q)

    global params
    mx   = params.mx;
    meqn = params.meqn;

    % this is the DO NOTHING explicit term
    f = zeros( mx, meqn );
    q = reshape( q, mx, meqn );
    f = ConstructL( q, t );
    f = reshape( f, mx*meqn, 1 );

end
