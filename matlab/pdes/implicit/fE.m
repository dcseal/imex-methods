% Explicit right hand side function
function f = fE(t,q)

    global params
    mx   = params.mx;
    meqn = params.meqn;

    % this is the DO NOTHING explicit term
    q = reshape( q, mx, meqn );
    f = ConstructL( q, t );
    f = zeros( size(q) );

end
