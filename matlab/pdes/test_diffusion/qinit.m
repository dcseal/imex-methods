% User supplied Initial condition function
function [q] = qinit(a,b)

    q = myquad( @q4func, a, b );

end

function q1 = q1func(x)

    q1 = sin(2*pi*x);

    lambda = 2*pi;
    q1 = sin( lambda*(x) );

end

function q2 = q2func(x)

    width = 0.2;
    I = find( x > 0.3 & x < 0.5 );
    q2 = zeros( size(x) );
    q2(I) = 1 - sin( pi*( x(I)-0.4 ) / (width) ) .^6 ;

end

function q1 = q3func(x)

    q1 = 2.0 + sin(4*x);

end

function q = q4func(x)

    global params;
    epsilon = params.eps;

    q = x.*(1-x).*exp(x);

tf = 0;
q = exp(-tf)*x.*(1-x);

end
