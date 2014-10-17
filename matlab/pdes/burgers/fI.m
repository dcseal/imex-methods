% Implicit right hand side function
function f = fI(t,q)

    global params
    mx   = params.mx;
    meqn = params.meqn;
    dx   = params.dx;
    mbc  = params.mbc;

    b = ones(mx, 1);

%   if( params.sorder == 4 )
%   
%       %%%% This section works for periodic and fourth order accuracy %%%%
%       Am = spdiags( [b -15*b 15*b -b], [-2; -1; 0; 1], mx, mx );
%       Ap = spdiags( [b -15*b 15*b -b], [-1; 0; 1; 2],  mx, mx );

%       params.AI = 1/(12*dx^2)*(Ap-Am);

%       row = params.AI(3, 1:5)';

%       % upper right block
%       params.AI(1, mx-1:mx) = row(1:2);
%       params.AI(2, mx)      = row(1);

%       % lower left block
%       params.AI(mx,   1:2)  = row(4:5);
%       params.AI(mx-1, 1)    = row(5);

%   elseif( params.sorder == 2 )

%       % Periodic and 2nd order accuracy
%       Am = spdiags( [ -b b ], [ -1; 0],  mx, mx );
%       Ap = spdiags( [ -b b ], [  0; 1],  mx, mx );

%       params.AI = 1/(dx^2)*(Ap-Am);

%       % upper left block
%       params.AI(1, mx) = ( 1  ) / dx^2;

%       % lower left block
%       params.AI(mx, 1) = ( 1  ) / dx^2;

%   end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Dirichlet and 2nd order Accuracy %%
%   Am = spdiags( [ -b b ], [ -1; 0],  mx, mx );
%   Ap = spdiags( [ -b b ], [  0; 1],  mx, mx );

%   params.AI = 1/(dx^2)*(Ap-Am);

%   % upper left block
%   params.AI(1, 1) = (-7/2) / dx^2;
%   params.AI(1, 2) = (-1/2) / dx^2;
%   params.AI(1, 3) = ( 1  ) / dx^2;

%   % lower left block
%   params.AI(mx, mx  ) = (-7/2) / dx^2;
%   params.AI(mx, mx-1) = (-1/2) / dx^2;
%   params.AI(mx, mx-2) = ( 1  ) / dx^2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Dirichlet and 4th order Accuracy %%
    Am = spdiags( [b -15*b 15*b -b], [-2; -1; 0; 1], mx, mx );
    Ap = spdiags( [b -15*b 15*b -b], [-1; 0; 1; 2],  mx, mx );

    params.AI = 1/(12*dx^2)*(Ap-Am);

    % upper left block
    params.AI(1, 1) =  481 * (-1/(72*dx^2)); 
    params.AI(1, 2) = -215 * (-1/(72*dx^2));
    params.AI(1, 3) =   37 * (-1/(72*dx^2));
    params.AI(1, 4) =   -3 * (-1/(72*dx^2));

    params.AI(2, 1) =  269 / (144*dx^2);
    params.AI(2, 2) = -403 / (144*dx^2);
    params.AI(2, 3) =  209 / (144*dx^2);
    params.AI(2, 4) =  -15 / (144*dx^2);

    % lower left block
    params.AI(mx, mx)   =  481 * (-1/(72*dx^2)); 
    params.AI(mx, mx-1) = -215 * (-1/(72*dx^2));
    params.AI(mx, mx-2) =   37 * (-1/(72*dx^2));
    params.AI(mx, mx-3) =   -3 * (-1/(72*dx^2));

    params.AI(mx-1, mx)   =  269 / (144*dx^2);
    params.AI(mx-1, mx-1) = -403 / (144*dx^2);
    params.AI(mx-1, mx-2) =  209 / (144*dx^2);
    params.AI(mx-1, mx-3) =  -15 / (144*dx^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    f = params.eps*params.AI*q;
    f = reshape(f, mx*meqn, 1 );

end
