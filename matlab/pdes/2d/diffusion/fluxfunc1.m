% Flux function for the problem
function F = fluxfunc1( q )

    global params

    F = params.u * q;

end
