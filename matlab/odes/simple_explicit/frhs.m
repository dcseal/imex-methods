% right hand side function (used for exponential integrator)
function f = frhs(q)

    global params

    % source term gets treated implicitly.
    f = params.lambda * q;

end
