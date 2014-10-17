% Explicit right hand side function
function f = fE(t,q)

    global params

    f    = zeros( size(q) );
    f(1) = -q(2);
    f(2) = q(1);

end
