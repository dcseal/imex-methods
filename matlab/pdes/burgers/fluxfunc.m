% Flux function for the problem
function F = fluxfunc( x, q )

    global params

    F = 0.5 * q.*q;

end
