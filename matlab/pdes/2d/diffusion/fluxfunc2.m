% Flux function for the problem
function F = fluxfunc2( x, q )

    global params

    F = params.v * q;

end
