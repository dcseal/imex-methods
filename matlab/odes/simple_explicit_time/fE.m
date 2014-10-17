% Explicit right hand side function
function f = fE(t,q)

    global params

    % source term gets treated implicitly.
    f = zeros( size(q) );
    f(1) = sin(t);

    % second variable represents time
    f(2) = 1;

end
