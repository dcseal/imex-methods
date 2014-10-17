% Flux function for the problem
function F = fluxfunc( x, q )

    global params

    F = params.u * q;

end
