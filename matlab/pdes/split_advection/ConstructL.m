function L = ConstructL( q, t )

    global params
    mx     = params.mx;
    node   = params.node;
    sorder = params.sorder;
    dx     = params.dx;
    mbc    = params.mbc;
    meqn   = params.meqn;


    [ql, qr] = set_boundary( q, t );
    qbig     = [ql; q; qr];

    % compute the function values at each interface
    I = 1:(mx+1);    
    F = zeros( mx+1, meqn );


    if( sorder == 2 )

        % two point symmetric stencil
        qim1 = qbig( I  , : );
        qi   = qbig( I+1, : );

        qimh = 0.5*( qim1 + qi );

    elseif( sorder == 3 )

        % TODO - THIS DEPENDS ON THE SPEED, SO THIS MAY BE THE INCORRECT FLUX
        % FUNCTION!

        % three point 'upwind' stencil.
        qim2 = qbig( I  , : );
        qim1 = qbig( I+1, : );
        qi   = qbig( I+2, : );

        qimh = 1/6*( 2*qi + 5*qim1 - qim2 );

    elseif (sorder == 4 )


        % four point symmetric stencil
        qim2 = qbig( I  , : );
        qim1 = qbig( I+1, : );
        qi   = qbig( I+2, : );
        qip1 = qbig( I+3, : );

        qimh = 7/12*( qi + qim1 ) - 1/12 * (qim2 + qip1);

    end


    %%%%%%%%%%%%%%%% At some point this worked %%%%%%%%%%%%   
%   for i=1:(mx+1)

%       % node point ( x_{i-1/2} )

%       im1 = i-1+mbc;
%       im2 = i-2+mbc;
%       ip1 = i+1+mbc;

%       if( sorder == 2 )

%           qi = 0.5*( qbig(iq,:) + qbig(im1,:) );

%       elseif( sorder == 3 )

%           qi = 1/6*( 2*qbig(iq,:) + 5*qbig(im1,:) - qbig(im2,:) );

%       elseif (sorder == 4 )

%           qi = 7/12*( qbig(iq,:) + qbig(im1,:) ) - 1/12*(qbig(im2,:) \
%                   + qbig(ip1,:) );

%       end

%       F(i,:) = fluxfunc( node(i), qi );

%   end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    F    = fluxfunc( node, qimh );

    Fp = F(2:end,:);
    Fm = F(1:mx, :);

    L = -1/dx*( Fp - Fm );

end
