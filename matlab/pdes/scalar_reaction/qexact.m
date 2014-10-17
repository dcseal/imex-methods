% User supplied Initial condition function
% function q = qexact(a,b)
function q = qexact(x)

    global params
    mx     = params.mx;
    meqn   = params.meqn;

    q = q3func(x);

%   q = quad( @q3func, a, b );

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

     flive = -4*pi*tau*cos(pi*x).^2/(1+4*pi^2*tau^2)+2*tau*pi/(1+4*pi^2*tau^2)+2*sin(pi*x).*cos(pi*x)/(1+4*pi^2*tau^2);

     fdecay = ...
     2*(4*pi*tau*cos(pi*x).^2.*cos(pi*t)^2-2*pi*tau*cos(pi*x).^2-2*pi*tau*cos(pi*t).^2+pi*tau+4*pi*tau*sin(pi*x).*cos(pi*x).*sin(pi*t).*cos(pi*t)-2*sin(pi*x).*cos(pi*x).*cos(pi*t).^2+sin(pi*x).*cos(pi*x)+2*cos(pi*x).^2*sin(pi*t).*cos(pi*t)-sin(pi*t).*cos(pi*t)).*exp(-t/tau)/(1+4*pi^2*tau^2);


%    fdecay = ...
%    2.*(4.*pi.*tau.*cos(pi.*x).^2.*cos(pi.*t).^2-2.*pi.*tau.*cos(pi.*x).^2-2.*pi.*tau.*cos(pi.*t).^2+pi.*tau+4.*pi.*tau.*sin(pi.*x).*cos(pi.*x).*sin(pi.*t).*cos(pi.*t)-2.*sin(pi.*x).*cos(pi.*x).*cos(pi.*t).^2+sin(pi.*x).*cos(pi.*x)+2.*cos(pi.*x).^2.*sin(pi.*t).*cos(pi.*t)-sin(pi.*t).*cos(pi.*t)).*exp(-t/tau)/(1+4.*pi.^2.*tau.^2);

    % This needs to be consistent with qinit %
    x = x - t;   % this really should be x-ut

    width = 0.6;
%   I     = find( x > 0.3 & x < 0.9 );
%   q2    = zeros( size(x) );
%   q2(I) = cos( pi*( x(I)-0.6 ) / (width) ) .^6 ;
    q2 = cos( 2*pi*(x) );

    q3 = flive + fdecay + exp(-t/tau) * q2;

end
