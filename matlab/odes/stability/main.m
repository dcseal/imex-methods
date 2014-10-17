%%% Main program.  See set_params for user defined options. %%%

    clear all;

    addpath('../lib/');

    mx = 100;
    my = 100;

    xmin = -4.5;
    xmax =  2.5;
    ymin = -5;
    ymax =  5;

    x = linspace(xmin, xmax, mx );
    y = linspace(ymin, ymax, my );
    [XX YY] = meshgrid( x, y );

    global params

    params.coeffs = 'RK3';

    ZZ = XX + sqrt(-1) * YY;

    params.mx = mx;
    params.my = my;
    params.meqn = mx*my;

    params.zz = zeros( mx*my, 1 );
    for i=1:mx
    for j=1:my
        params.zz( (i-1)*mx + j ) = ZZ(i,j)
    end
    end


    for j=2:5

        figure(j-1); clf;

        if( j == 2 )
            params.coeffs = 'RK2';
        elseif( j == 3 )
            params.coeffs = 'RK3';
        elseif( j == 4 )
            params.coeffs = 'RK4';
        else
            params.coeffs = 'dumb';
        end

        % Take a single time step of length dt = 1.0 with the integrator
        tn   = 0; dt = 1.0;

        params.sdc_order = j;
        qnp1 = sdc_integrator( tn, dt, ones(1, mx*my) );

        qbig = zeros( mx, my );
        for j=1:my
        for i=1:mx
            qbig( i,j ) = qnp1( (i-1)*mx + j );
        end
        end
        qnp1 = qbig;

        hold on;
        if( j == 2 )
            contour( XX, YY, 1 - abs( qnp1 ), [0 0], '-g', 'LineWidth', 3 );
        elseif( j == 3 )
            contour( XX, YY, 1 - abs( qnp1 ), [0 0], '-r', 'LineWidth', 3 );
        else
            contour( XX, YY, 1 - abs( qnp1 ), [0 0], '-b', 'LineWidth', 3 );
        end
            
        hold on;
        plot( [xmin xmax], [0 0], '-k' );
        hold on;
        plot( [0 0 ], [ymin ymax], '-k' );
        axis( [xmin xmax ymin ymax] );
        hold off;

        t = title('Stability Regions for Explicit Forward Euler SDC Methods');

    end

%   l = legend('2nd Order', '3rd Order', '4th Order', 'Location', 'NorthWest' );
%   hold off;

