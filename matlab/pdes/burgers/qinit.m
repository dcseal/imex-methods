% User supplied Initial condition function
function [q] = qinit(a,b)

    global params
    mx     = params.mx;
    meqn   = params.meqn;

    q = quad( @q1func, a, b );

end

function q1 = q1func(x)

    global params
    q1 = ( 2*pi*params.eps*sin(pi*x) ) ./ ( params.a + cos(pi*x) );

end
