% Explicit right hand side function
function f = fE(t,q)

    global params
    f = zeros( size( q ) );
    f = reshape( q'.*params.zz, 1, params.meqn );

end
