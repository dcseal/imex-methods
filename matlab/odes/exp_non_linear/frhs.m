% right hand side function (used for exponential integrator)
function f = frhs(q)

    global params

    f = zeros( size(q) );
    f(1) = -q(1)^2;
    f(2) = 1.0;


end
