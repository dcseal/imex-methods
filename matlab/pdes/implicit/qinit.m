% User supplied Initial condition function
function [q] = qinit(a,b)

    global params
    mx     = params.mx;
    meqn   = params.meqn;

    q = quad( @q3func, a, b );

end

function q1 = q1func(x)

    q1 = sin(2*pi*x);

    lambda = 2*pi;
    q1 = sin( lambda*(x) ) + 0.01*cos( lambda*(x) );

end

function q2 = q2func(x)

    width = 0.2;
    I = find( x > 0.3 & x < 0.5 );
    q2 = zeros( size(x) );
    q2(I) = 1 - sin( pi*( x(I)-0.4 ) / (width) ) .^6 ;

end

function q1 = q3func(x)

    q1 = ones( size(x) );

end

