% User supplied Initial condition function
function q = qexact(a,b)

    global params
    mx     = params.mx;
    meqn   = params.meqn;

    q = quad( @qfunc, a, b );

end

function q = qfunc(x)

    width = 0.2;
    I = find( x > 0.5 & x < 0.7 );
    q = zeros( size(x) );
    q(I) = 1 - sin( pi*( x(I)-0.4 ) / (width) ) .^6 ;

end


