function L = ConstructL( q, t )

    global params
    mx     = params.mx;
    my     = params.my;
    node   = params.node;
    sorder = params.sorder;
    dx     = params.dx;
    dy     = params.dy;
    mbc    = params.mbc;
    meqn   = params.meqn;

    qbig     = fill_boundary( q, t);

    % compute the function values at each interface
    I  = 1:(mx);    Ic = I+mbc;
    J  = 1:(my);    Jc = J+mbc;
    Il = 1:(mx+1);    
    Jl = 1:(my+1);    
    if( sorder == 2 )

        % two point symmetric stencil
        qim1 = qbig( Il  , J+mbc, : );
        qi   = qbig( Il+1, J+mbc, : );

        % two point symmetric stencil
        qjm1 = qbig( I+mbc,  Jl,   : );
        qj   = qbig( I+mbc,  Jl+1, : );

        % face averages for stencils in the x and y direction
        qimh = 0.5*( qim1 + qi );
        qjmh = 0.5*( qjm1 + qj );


%   elseif( sorder == 3 )

%       % TODO - THIS DEPENDS ON THE SPEED, SO THIS MAY BE THE INCORRECT FLUX
%       % FUNCTION!

%       % three point 'upwind' stencil.
%       qim2 = qbig( I  , : );
%       qim1 = qbig( I+1, : );
%       qi   = qbig( I+2, : );

%       qimh = 1/6*( 2*qi + 5*qim1 - qim2 );

    elseif (sorder == 4 )


        % four point symmetric stencil
        qim2 = qbig( Il  , Jc, : );
        qim1 = qbig( Il+1, Jc, : );
        qi   = qbig( Il+2, Jc, : );
        qip1 = qbig( Il+3, Jc, : );

        qimh = 7/12*( qi + qim1 ) - 1/12 * (qim2 + qip1);

        % four point symmetric stencil
        qjm2 = qbig( Ic, Jl  , : );
        qjm1 = qbig( Ic, Jl+1, : );
        qj   = qbig( Ic, Jl+2, : );
        qjp1 = qbig( Ic, Jl+3, : );

        qjmh = 7/12*( qj + qjm1 ) - 1/12 * (qjm2 + qjp1);


    end

    F  = fluxfunc1( qimh );
    Fp = F(2:end, J, :);
    Fm = F(1:mx,  J, :);

    G  = fluxfunc2( node, qjmh );
    Gp = G(I, 2:end, :);
    Gm = G(I, 1:my,  :);

    L = -1/dx*( Fp - Fm ) - 1/dy*( Gp - Gm );

end
