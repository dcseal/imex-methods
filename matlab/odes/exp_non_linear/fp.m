% Derivative of rhs function.  This is used for the exponential integrator
function A = fp(q)

    global params;

    A = [ -2*q(1) 0; 0 0 ];

end
