% Explicit right hand side function
function f = fE(t,q)

    global params

    % source term gets treated implicitly.
    f = zeros( size(q) );
    f(1) = -2*q(1);
    f(2) = q(1) - q(2)*(1+q(2));

end
