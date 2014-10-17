% User supplied Initial condition function
function q = qexact(a,b)

    global params
    mx     = params.mx;
    meqn   = params.meqn;

    q = myquad( @q3func, a, b, 1e-12 );

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

function q3 = q3func(x)

    global params
    tau = params.tau;
    t   = params.tf;

    fdecay = ...
    -2.*(tau.^2.*sin(pi.*t).*cos(pi.*t)+8.*tau.^3.*pi.^3.*cos(pi.*t).^2+2.*tau.*pi.*cos(pi.*t).^2-sin(pi.*x).*cos(pi.*x)+2.*sin(pi.*x).*cos(pi.*x).*cos(pi.*t).^2-2.*tau.^3.*pi.*cos(pi.*t).^2-2.*cos(pi.*x).^2.*sin(pi.*t).*cos(pi.*t)-8.*tau.^2.*pi.^2.*cos(pi.*x).^2.*sin(pi.*t).*cos(pi.*t)-4.*tau.*pi.*sin(pi.*x).*cos(pi.*x).*sin(pi.*t).*cos(pi.*t)+8.*tau.^2.*pi.^2.*sin(pi.*x).*cos(pi.*x).*cos(pi.*t).^2-16.*tau.^3.*pi.^3.*sin(pi.*x).*cos(pi.*x).*sin(pi.*t).*cos(pi.*t)+4.*tau.^3.*pi.*sin(pi.*x).*cos(pi.*x).*sin(pi.*t).*cos(pi.*t)+sin(pi.*t).*cos(pi.*t)-2.*tau.^2.*cos(pi.*x).^2.*sin(pi.*t).*cos(pi.*t)+4.*tau.^3.*pi.*cos(pi.*x).^2.*cos(pi.*t).^2-16.*tau.^3.*pi.^3.*cos(pi.*x).^2.*cos(pi.*t).^2-4.*tau.*pi.*cos(pi.*x).^2.*cos(pi.*t).^2+4.*tau.^2.*pi.^2.*sin(pi.*t).*cos(pi.*t)+2.*tau.^2.*sin(pi.*x).*cos(pi.*x).*cos(pi.*t).^2-4.*tau.^2.*pi.^2.*sin(pi.*x).*cos(pi.*x)-tau.*pi+2.*tau.*pi.*cos(pi.*x).^2-4.*tau.^3.*pi.^3+8.*tau.^3.*pi.^3.*cos(pi.*x).^2-tau.^2.*sin(pi.*x).*cos(pi.*x)+tau.^3.*pi-2.*tau.^3.*pi.*cos(pi.*x).^2)/(1+8.*tau.^2.*pi.^2+2.*tau.^2+16.*tau.^4.*pi.^4-8.*tau.^4.*pi.^2+tau.^4);

    flong = ...
    -2.*(-4.*tau.^2.*pi.^2.*cos(t).*sin(pi.*x).*cos(pi.*x)-tau.^2.*cos(t).*sin(pi.*x).*cos(pi.*x)+tau.^3.*pi.*cos(t)-2.*tau.^2.*pi.*sin(t)-tau.^3.*sin(t).*sin(pi.*x).*cos(pi.*x)+8.*tau.^3.*pi.^3.*cos(t).*cos(pi.*x).^2-2.*tau.^3.*pi.*cos(t).*cos(pi.*x).^2+2.*tau.*pi.*cos(t).*cos(pi.*x).^2-tau.*sin(t).*sin(pi.*x).*cos(pi.*x)+4.*tau.^2.*pi.*sin(t).*cos(pi.*x).^2-4.*tau.^3.*pi.^3.*cos(t)+4.*tau.^3.*pi.^2.*sin(t).*sin(pi.*x).*cos(pi.*x)-tau.*pi.*cos(t)-cos(t).*sin(pi.*x).*cos(pi.*x))/(1+8.*tau.^2.*pi.^2+2.*tau.^2+16.*tau.^4.*pi.^4-8.*tau.^4.*pi.^2+tau.^4);

    % This needs to be consistent with qinit %
    width = 0.6;
    x = x - t;
    I     = find( x > 0.3 & x < 0.9 );
    q2    = zeros( size(x) );
    q2(I) = cos( pi*( x(I)-0.6 ) / (width) ) .^6 ;

    q3 = flong + exp(-t/tau)*( fdecay + q2 );

end
