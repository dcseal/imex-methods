% User supplied Initial condition function
function q = qexact(xpts)

    global params
    mx     = params.mx;
    meqn   = params.meqn;

    lambda = 2*pi;
    tf = params.tfinal;

    [mpts, mdim] = size(xpts);
    q = zeros( mpts, 1 );
    for n=1:mpts

        x = xpts(n,1);
        y = xpts(n,2);

        q(n,1) = sin(lambda*y) * sin( lambda*(x- params.u*tf) ) .* exp( - params.eps * lambda^2 * tf );

    end

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
