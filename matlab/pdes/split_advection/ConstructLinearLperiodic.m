% Periodic Boundary Conditions
%
% This function constructs the matrix A such that in the ODE
% q' = -1/dx( f(q_{i+1/2} - f( q_{i-1/2} ) ) ), the right hand side is given by
% q' = A q.
%
function A = ConstructLinearLperiodic( params )

    mx     = params.mx;
    dx     = params.dx;
    meqn   = params.meqn;
    sorder = params.sorder;
    u      = params.u;

    if( sorder == 2 )
        
        b = ones( mx, meqn );
        A = spdiags( [-b b], [-1; 1], mx, mx );
        
        % set boundary values:
        A(meqn*mx-(1:meqn), 1 ) = 1;
        A(1,mx) = -1;

        A = -u/(2*dx)*A;

    end

    if( sorder == 3 )
        
        b = ones(mx,1);
        A = spdiags( [b -6*b 3*b 2*b], [-2; -1; 0; 1], mx, mx );
        
        % first row
        A(1,mx-1) = 1;
        A(1,mx)   = -6;

        % second row
        A(2,mx) = 1;

        % last row
        A(mx,1) = 2;

        A = -u/(6*dx)*A;

    end

    if( sorder == 4 )
        
        b = ones(mx,1);
        A = spdiags( [b -8*b 8*b -b], [-2; -1; 1; 2], mx, mx );
        
        % first row
        A(1,mx-1) = 1;
        A(1,mx)   = -8;

        % second row
        A(2,mx) = 1;

        % second to last row:
        A(mx-1, 1) = -1;

        % last row
        A(mx, 1) = 8;
        A(mx, 2) = -1;

        A = -u/(12*dx)*A;

    end




end
