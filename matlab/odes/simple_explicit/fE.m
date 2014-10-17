% Explicit right hand side function
function f = fE(t,q)

    global params

    % source term gets treated implicitly.
    f = params.lambda * q;

end
