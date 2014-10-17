% User supplied Initial condition function
% function [q] = qinit(a,b)
function [q] = qinit(x)

    global params
    mx     = params.mx;
    meqn   = params.meqn;

    q = q1func(x);
    %q = quad( @q1func, a, b );

end

function q1 = q1func(x)

    q1 = cos(2*pi*x);

end

function q2 = q2func(x)

    width = 0.6;
    I     = find( x > 0.3 & x < 0.9 );
    q2    = zeros( size(x) );
    q2(I) = cos( pi*( x(I)-0.6 ) / (width) ) .^6 ;

end

function q = q3func(x)

    q = zeros(size(x));

end
