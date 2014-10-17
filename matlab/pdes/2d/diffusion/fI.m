% Implicit right hand side function
function f = fI(t,q)

    global params
    mx   = params.mx;
    my   = params.my;

    meqn = params.meqn;
    dx   = params.dx;
    dy   = params.dy;

    mbc  = params.mbc;

    b = ones(mx, 1);

    %%%% This section works for periodic and fourth order accuracy %%%%
    Am = spdiags( [b -15*b 15*b -b], [-2; -1; 0; 1], mx, mx );
    Ap = spdiags( [b -15*b 15*b -b], [-1; 0; 1; 2],  mx, mx );

    Axx = 1/(12*dx^2)*(Ap-Am);

    row = Axx(3, 1:5)';

    % upper right block
    Axx(1, mx-1:mx) = row(1:2);
    Axx(2, mx)        = row(1);

    % lower left block
    Axx(mx,   1:2)  = row(4:5);
    Axx(mx-1, 1)    = row(5);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Dirichlet and 2nd order Accuracy %%
%   Am = spdiags( [ -b b ], [ -1; 0],  mx, mx );
%   Ap = spdiags( [ -b b ], [  0; 1],  mx, mx );

%   Axx = 1/(dx^2)*(Ap-Am);

%   % upper left block
%   Axx(1, 1) = (-7/2) / dx^2;
%   Axx(1, 2) = (-1/2) / dx^2;
%   Axx(1, 3) = ( 1  ) / dx^2;

%   % lower left block
%   Axx(mx, mx  ) = (-7/2) / dx^2;
%   Axx(mx, mx-1) = (-1/2) / dx^2;
%   Axx(mx, mx-2) = ( 1  ) / dx^2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Dirichlet and 4th order Accuracy %%
%   Am = spdiags( [b -15*b 15*b -b], [-2; -1; 0; 1], mx, mx );
%   Ap = spdiags( [b -15*b 15*b -b], [-1; 0; 1; 2],  mx, mx );

%   Axx = 1/(12*dx^2)*(Ap-Am);

%   % upper left block
%   Axx(1, 1) =  481 * (-1/(72*dx^2)); 
%   Axx(1, 2) = -215 * (-1/(72*dx^2));
%   Axx(1, 3) =   37 * (-1/(72*dx^2));
%   Axx(1, 4) =   -3 * (-1/(72*dx^2));

%   Axx(2, 1) =  269 / (144*dx^2);
%   Axx(2, 2) = -403 / (144*dx^2);
%   Axx(2, 3) =  209 / (144*dx^2);
%   Axx(2, 4) =  -15 / (144*dx^2);

%   % lower left block
%   Axx(mx, mx)   =  481 * (-1/(72*dx^2)); 
%   Axx(mx, mx-1) = -215 * (-1/(72*dx^2));
%   Axx(mx, mx-2) =   37 * (-1/(72*dx^2));
%   Axx(mx, mx-3) =   -3 * (-1/(72*dx^2));

%   Axx(mx-1, mx)   =  269 / (144*dx^2);
%   Axx(mx-1, mx-1) = -403 / (144*dx^2);
%   Axx(mx-1, mx-2) =  209 / (144*dx^2);
%   Axx(mx-1, mx-3) =  -15 / (144*dx^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % place the matrix along the main diagonal
    B = speye( my );
    params.AI = kron( B, Axx );


    f = params.eps*params.AI*q;
    f = reshape(f, mx*my*meqn, 1 );

end
