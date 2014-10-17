% Flux function for the problem
function F = fluxfunc( xi, qi )

    global params
    F = params.u * qi;

end
