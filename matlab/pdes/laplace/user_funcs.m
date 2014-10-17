function [qex, q, aux] = user_funcs(a,b)

    global params
    mx     = params.mx;
    meqn   = params.meqn;

    lambda = 2*pi;

    qex = myquad( @q1func, a, b );
    q   = myquad( @q1func0, a, b );
    aux = myquad( @auxinit, a, b );


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


function q1 = q1func0(x)

    q1 = sin(2*pi*x);

    lambda = 2*pi;
    q1 = sin( lambda*(x) );

end

% User supplied Initial condition function
function a = auxinit(x)

    width = 0.2;
    I = find( x > 0.3 & x < 0.5 );
    a = zeros( size(x) );
    a(I) = 1 - sin( pi*( x(I)-0.4 ) / (width) ) .^6 ;

end
