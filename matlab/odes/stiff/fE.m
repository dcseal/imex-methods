% Explicit right hand side function
function f = fE(t,y)

    global params;

    % source term gets treated implicitly.
    f = cos(t) / params.tau;

end
