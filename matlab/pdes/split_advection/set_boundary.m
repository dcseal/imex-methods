function [ql, qr] = set_boundary( q, t )

    global params

    mbc    = params.mbc;
    meqn   = params.meqn;
    sorder = params.sorder;
    mx     = params.mx;

    ql = zeros(mbc,meqn);
    qr = zeros(mbc,meqn);
    if( strcmp( params.mbc_type_left , 'periodic' ) )

        % Trivial copy
        for i=1:mbc

            ql(i, :)  = q(mx-mbc+i, :);
            qr(i, :)  = q(i, : );

        end

    end

    % dirichlet/inflow/outflow boundary conditions
    if( strcmp( params.mbc_type_left, 'dirichlet' ) )

        if( t > 0.5 & t < 0.7 )
            width = 0.2;
            A = cos( pi*( t-0.4 ) / (width) ) .^6 ;
        else
            A = 0;
        end

        if( mbc == 1 )

            % inflow
            ql(1,:) = 2*A-q(1,:);

            % outflow
            qr(1,:) = 2*q(mx,:) - q(mx-1,:);

        elseif( mbc == 2 )

            % inflow
            ql(2,:) =  4*A - 13/3*q(1,:) +  5/3*q(2,:)-1/3*q(3,:);
            ql(1,:) = -12*A + 7*(q(1,:)+ql(2,:)) - q(2,:);

        end

    elseif( strcmp( params.mbc_type_left, 'outflow' ) )

        disp('you need to be written!');

    end

    if( strcmp( params.mbc_type_right, 'dirichlet' ) )

        B = 0;

        % diriclet!
        qr(1,:) = 4*B +1/3*(-13*q(mx,:) + 5*q(mx-1,:) -   q(mx-2,:) );
        qr(2,:) = 16*B+1/3*(-70*q(mx,:) +32*q(mx-1,:) - 7*q(mx-2,:) );


    elseif( strcmp( params.mbc_type_right, 'outflow' ) )

        % outflow
        % conditions from my maple script
        qr(1,:) = -q(mx-3,:)+4*q(mx-2,:)-6*q(mx-1,:)+4*q(mx,:);
        qr(2,:) = -q(mx-2,:)+4*q(mx-1,:)-6*q(mx,:)+4*qr(1,:);


        % stencil borrowed from:
        % A high-order finite-volume method for hyperbolic conservation laws
        % on locally-refined grids
        % ( https://seesar.lbl.gov/anag/publications/petermc/HOFV2010.pdf )
        %           qr(1,:) =(-3*q(mx-3,:) +13*q(mx-2,:)- ...
                %                     23*q(mx-1,:) +25*q(mx,:) )/12;

        %           qr(2,:) =(  q(mx-2,:)- 5*q(mx-1,:)+ ...
                %                     13*q(mx,:)+3*qr(1,:) )/12;



        %           THIS IS BROKEN or perhaps not stable
        %           qr(1,:) =  3*q(mx-3,:) -16*q(mx-2,:)+ ...
        %                     22*q(mx-1,:) - 8*q(mx,:);

        %           qr(2,:) = 3*q(mx-2,:)-16*q(mx-1,:)+ ...
        %                     22*q(mx,:)-8*qr(1,:);



    end

end
