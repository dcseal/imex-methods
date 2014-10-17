% User supplied Initial condition function
% function q = qexact(a,b)
function q = qexact(x)

    global params
    mx     = params.mx;
    meqn   = params.meqn;

    q = zeros( length(x), 2 );
    q(:,1) = q1func(x);
    q(:,2) = q2func(x);

end

function q1 = q1func(x)

    q1 = cos(2*pi*x);

end

function q2 = q2func(x)

    width = 0.2;
    I = find( x > 0.3 & x < 0.5 );
    q2 = zeros( size(x) );
    q2(I) = cos( pi*( x(I)-0.4 ) / (width) ) .^6 ;

end

