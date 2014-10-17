% right hand side function (used for exponential integrator)
function f = frhs(q)

    f    = zeros( size(q) );
    f(1) = -q(2);
    f(2) = q(1);


end
