function L = ConstructL( q, t )
%CONSTRUCTL.   Construct the right hand side for q_t = L( q ).
%
% The routine L = CONSTRUCTL( q, t ) constructs the RHS function for a finite volume method.  It is
% assumed that the right hand side is of the conservative flux difference form
%
%       q_t = L( q ) = -1/dx ( f( q_{i+1/2} ) - f( q_{i-1/2} ) )
%
% for the PDE q_t + f(q)_x = 0.
%
% The shape of L is the same as of q, so that an ODE solver can be applied to
% the result.
%
% Required routines (that this function should know exist:
%
%       set_boundary - sets the boundary conditions on the solution
%       fluxfunc     - the flux function
%
% Input
% -----
%
%       q    - current solution (total of mx elements)
%       t    - current time ( scalar )
%
% Returns
% -------
%
%       L  - Right hand side function q_t = L( q ).

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

    elseif (sorder == 6 )

        % four point symmetric stencil
        qim3 = qbig( I  , : );
        qim2 = qbig( I+1, : );
        qim1 = qbig( I+2, : );
        qi   = qbig( I+3, : );
        qip1 = qbig( I+4, : );
        qip2 = qbig( I+5, : );

        qimh = 1/60*( qim3 + qip2 ) + 37/60 * (qi + qim1 ) - 2/15 * ( qim2 + qip1 );

    end

    F  = fluxfunc( node, qimh );

    Fp = F(2:end,:);
    Fm = F(1:mx, :);

    L = -1/dx*( Fp - Fm );

end
