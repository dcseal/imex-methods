% Inflow/Outflow Boundary Conditions
%
% This routine constructs the matrix A with inflow boundary conditions on the
% left hand side, and outflow boundary conditions on the right hand side.  That
% is, q_{-1/2} = A, where A is a user defined function.
%
%function A = ConstructLinearLinflow( mx, dx, meqn, sorder, u)
function A = ConstructLinearLinflow( param )

    mx     = param.mx;
    dx     = param.dx;
    meqn   = param.meqn;
    sorder = param.sorder;
    u      = param.u;

    if( sorder == 2 )
        
        b = ones(mx,1);
        A = spdiags( [-b b], [-1; 1], mx, mx );

        % Inflow part of the boundary:
        %        
        % note: to second order, we can say that q_0 = 2*A-q_1.
        % This implies that
        %
        %       \partial_t q_1 = -u/(2dx) (-q1 + q2 ) - A*u / dx.
        %
        A(1,1) = 1;
        A(1,2) = 1;

        % Outflow part of the boundary:
        %
        % Again, we can say that q_{mx+1/2} = 2 q(mx) - q(mx-1) 
        A(mx,mx-1) = -3;
        A(mx,mx)   = 3;

        A = -u/(2*dx)*A;

    end

    if( sorder == 3 )
        
%       b = ones(mx,1);
%       A = spdiags( [b -6*b 3*b 2*b], [-2; -1; 0; 1], mx, mx );
%       
%       % first row
%       A(1,mx-1) = 1;
%       A(1,mx)   = -6;

%       % second row
%       A(2,mx) = 1;

%       % last row
%       A(mx,1) = 2;

%       A = -u/(6*dx)*A;

    end

    if( sorder == 4 )
        
        b = ones(mx,1);
        A = spdiags( [b -8*b 8*b -b], [-2; -1; 1; 2], mx, mx );
       
        % inflow part of the boundary:

        % TODO

        % outflow part of the boundary:

        % TODO  

        % first row
%       A(1,mx-1) = 1;
%       A(1,mx)   = -8;

%       % second row
%       A(2,mx) = 1;

%       % second to last row:
%       A(mx-1, 1) = -1;

%       % last row
%       A(mx, 1) = 8;
%       A(mx, 2) = -1;

        A = -u/(12*dx)*A;

    end




end
