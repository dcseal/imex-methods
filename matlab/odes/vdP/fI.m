% Implicit right hand side function
function f = fI(t,q)

    global params
    eps = params.eps;

    f = zeros( size(q) );
    f(1) = 0.0;
    f(2) = 1/eps*( (1-q(1)^2) * q(2) - q(1) );

end
