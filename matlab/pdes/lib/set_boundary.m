function [ql, qr] = set_boundary( q, t )

    global params

    mbc    = params.mbc;
    meqn   = params.meqn;
    sorder = params.sorder;
    mx     = params.mx;

    ql = zeros(mbc,meqn);
    qr = zeros(mbc,meqn);

    % periodic means periodic!
    if( strcmp( params.mbc_type_left , 'periodic' ) )

        % Trivial copy
        for i=1:mbc

            ql(i, :)  = q(mx-mbc+i, :);
            qr(i, :)  = q(i, : );

        end

        return;

    end

    % boundary condition for hte left hand side
    if( strcmp( params.mbc_type_left, 'dirichlet' ) )

        A = 0.0;

        if( mbc == 1 )

            % inflow
            ql(1,:) = 2*A-q(1,:);

            % outflow
            qr(1,:) = 2*q(mx,:) - q(mx-1,:);

        elseif( mbc == 2 )

            % inflow
            ql(2,:) =  4*A - 13/3*q(1,:) +  5/3*q(2,:)-1/3*q(3,:);

            ql(1,:) = 16*A - 70/3*q(1,:) + 32/3*q(2,:)-7/3*q(3,:);
            %           ql(1,:) =-12*A + 7*( ql(2,:)+q(1,:) ) - q(2,:);

        end

    elseif( strcmp( params.mbc_type_left, 'outflow' ) )

        disp('you need to be written!');

    else
        disp('unrecognized mbc for left hand side');
    end

    % Boundary condition for the right hand side
    if( strcmp( params.mbc_type_right, 'dirichlet' ) )

        B = 0.0;

        % diriclet!
        qr(1,:) = 4*B +1/3*(-13*q(mx,:) + 5*q(mx-1,:) -   q(mx-2,:) );
        qr(2,:) = 16*B+1/3*(-70*q(mx,:) +32*q(mx-1,:) - 7*q(mx-2,:) );


    elseif( strcmp( params.mbc_type_right, 'outflow' ) )

        % outflow
        % conditions from my maple script
        qr(1,:) = -q(mx-3,:)+4*q(mx-2,:)-6*q(mx-1,:)+4*q(mx,:);
        qr(2,:) = -q(mx-2,:)+4*q(mx-1,:)-6*q(mx,:)+4*qr(1,:);

    else
        disp('unrecognized mbc for left hand side');
    end

end
