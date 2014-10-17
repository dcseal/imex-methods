% Implicit right hand side function
function f = fI(t,q)

    global params

    f = zeros( size(q) );
    f(1) = 1/params.eps*( q(2)^2 - q(1) );
    f(2) = 0;

end
