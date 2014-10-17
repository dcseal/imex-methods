function ynp1 = exp_integrator( tn, dt, yn )

    h  = dt;
    y0 = yn';

    % evaluate derivative of rhs
    A0 = fp( yn );
    f0 = frhs( y0' )';


    % method of Hochbruck et a. (Exp4)
    % M. Hochbruck, C. Lubich, H. Selhofer, ``Exponental Integrators for large
    % systems of differential equations'', SIAM Journal of Scientific Computing, 19,
    % (1998)

    % phi1(z) = ( exp(z) - 1 ) / z = 1 + z / 2! + z^2 / 3! + z^3 / 4! + z^4 / 5! + \cdots

    id = speye( size(A0) );

    % This should simply use one krylov solve

    nmax = 10;  % number of times to sum the small matrix when evaluating phi1( V) 
    n = 10;      % number of krylov vectors produced

    %%%%%%%%%% Compute A single Krylov Projection %%%%%%%%%%%%%
    [V,H] = arnoldi(h/3*A0,f0,n);

    % compute f( Hm ) = ( exp(Hm) - 1 ) / Hm 
    fm = eye( size(H) );
    for n=1:nmax
        fm = fm + H^n / factorial(n+1);
    end
    % compute ( exp( A0*h ) - 1 ) / (A0*h) * f0
    k1 = norm( f0, 2 ) * V * fm;
    k1 = k1(:,1);

    [V,H] = arnoldi(2*h/3*A0,f0,n);

    % compute f( Hm ) = ( exp(Hm) - 1 ) / Hm 
    fm = eye( size(H) );
    for n=1:10
        fm = fm + H^n / factorial(n+1);
    end
    % compute ( exp( A0*h ) - 1 ) / (A0*h) * f0
    k2 = norm( f0, 2 ) * V * fm;
    k2 = k2(:,1);

    [V,H] = arnoldi(h*A0,f0,n);

    % compute f( Hm ) = ( exp(Hm) - 1 ) / Hm 
    fm = eye( size(H) );
    for n=1:10
        fm = fm + H^n / factorial(n+1);
    end
    % compute ( exp( A0*h ) - 1 ) / (A0*h) * f0
    k3 = ( norm( f0, 2 ) * V * fm );
    k3 = k3(:,1);


    At = 1/3*h*A0;
    phi1 = id + At / 2 + At^2 / 6 + At^3 / 24 + At^4 / 120;


    At = 2/3*h*A0;
    phi2 = id + At / 2 + At^2 / 6 + At^3 / 24 + At^4 / 120;

    At = h*A0;
    phi3 = id + At / 2 + At^2 / 6 + At^3 / 24 + At^4 / 120;

    k1 = phi1 * f0;
    k2 = phi2 * f0;
    k3 = phi3 * f0;

    w4 = -7/300*k1 + 97/150*k2 - 37/300*k3;
    u4 = y0 + h*w4;
    r4 = frhs(u4')' - f0 - h*A0*w4;

    %  This should be one krylov solve here as well

    k4 = phi1*r4; k5 = phi2*r4; k6 = phi3*r4;

    w7 = 59/300*k1 - 7/75*k2 + 269/300*k3 + 2/3*( k4+k5+k6 );
    u7 = y0 + h*w7;
    r7 = frhs(u7')' - f0 - h*A0*w7;

    k7 = phi1*r7;

    ynp1 = y0 + h*(k3+k4-4/3*k5 + k6 + k7/6 );


    %%%%%%%%%%%%%% two-stage RK scheme (can be found in Hochbruck): %%%%%%%%%%%%%

    %%%%%%%%%% Compute A single Krylov Projection %%%%%%%%%%%%%
%   n = 2;  % number of krylov vectors produced
%   [V,H] = arnoldi(h*A0,f0,n);

%   % compute f( Hm ) = ( exp(Hm) - 1 ) / Hm 
%   fm = eye( size(H) );
%   for n=1:100
%       fm = fm + H^n / factorial(n+1);
%   end

%   % compute ( exp( A0*h ) - 1 ) / (A0*h) * f0
%   k1 = norm( f0, 2 ) * V * fm;
%   k1 = k1(:,1);

%   ynp1 = ( y0 + h*k1 )';

end
