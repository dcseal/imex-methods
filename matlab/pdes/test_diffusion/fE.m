% Explicit right hand side function
function f = fE(t,q)

    global params
    mx   = params.mx;
    meqn = params.meqn;
    aux  = params.aux;

    % advective plus source term
    f = zeros( mx, meqn );
    q = reshape( q, mx, meqn );
    % f = ConstructL( q, t) + exp(-t)*aux;
    f = exp(-t)*aux;
    f = reshape( f, mx*meqn, 1 );

end
