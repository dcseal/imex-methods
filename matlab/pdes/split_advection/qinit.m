% User supplied Initial condition function
function [q] = qinit(a,b)

    global params
    mx     = params.mx;
    meqn   = params.meqn;

    q = quad( @q1func, a, b );

end

function q1 = q1func(x)

    q1 = sin(2*pi*x);

end

function q2 = q2func(x)

    width = 0.2;
    I = find( x > 0.3 & x < 0.5 );
    q2 = zeros( size(x) );
    q2(I) = 1 - sin( pi*( x(I)-0.4 ) / (width) ) .^6 ;

%   q2 = zeros(size(x));
%   q2 = 0.5 + 0.25 * sin(2*pi*x);

end


