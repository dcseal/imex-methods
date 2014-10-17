% Explicit right hand side function
function f = fE(t,y)

    global params
    eps = params.eps;

    % source term gets treated implicitly.
    f = zeros( size(y) );
    f(1) = y(2);
    f(2) = 0.0;

end
