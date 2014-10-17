% User supplied Initial condition function
function q = qexact(a,b)

    global params
    mx     = params.mx;
    meqn   = params.meqn;

    q = quad( @q1func, a, b );

end

function q1 = q1func(x)
% q = sin( \lambda x )

    global params

    lambda = 2*pi;
    tf = params.tfinal;
    q1 = sin( lambda*(x- params.u*tf) ) .* exp( - params.eps * lambda^2 * tf );

end

function q2 = q2func(x)

    width = 0.2;
    I = find( x > 0.3 & x < 0.5 );
    q2 = zeros( size(x) );
    q2(I) = 1 - sin( pi*( x(I)-0.4 ) / (width) ) .^6 ;

end


function q = q3func(x)

    width = 0.2;
    I = find( x > 0.5 & x < 0.7 );
    q = zeros( size(x) );
    q(I) = 1 - sin( pi*( x(I)-0.4 ) / (width) ) .^6 ;

end
