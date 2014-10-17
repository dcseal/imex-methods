% Implicit right hand side function and construction of implicit matrix for the
% solver
function f = fI(t,q)

    global params
    mx   = params.mx;
    meqn = params.meqn;
    dx   = params.dx;
    mbc  = params.mbc;

    b = ones(mx, 1);


    if( strcmp( params.mbc_type_left, 'periodic' ) )


        if( mbc == 3 )

            %%% This section is supposed to provide 6th order accuracy %%%
            Am = spdiags( [-2*b 25*b -245*b 245*b -25*b 2*b], [-3; -2; -1; 0; 1; 2], mx, mx );
            Ap = spdiags( [-2*b 25*b -245*b 245*b -25*b 2*b], [-2; -1;  0; 1; 2; 3], mx, mx );
            params.AI = 1/(180*dx^2)*(Ap-Am);


            row = params.AI(4, 1:7)';

            % upper right block
            params.AI(1, mx-2:mx) = row(1:3);
            params.AI(2, mx-1:mx) = row(1:2);
            params.AI(3, mx)      = row(1);

            % lower left block
            params.AI(mx,   1:3)  = row(5:7);
            params.AI(mx-1, 1:2)  = row(6:7);
            params.AI(mx-2, 1)    = row(7);

        else

            %%%% This section works for periodic and fourth order accuracy %%%%
            Am = spdiags( [b -15*b 15*b -b], [-2; -1; 0; 1], mx, mx );
            Ap = spdiags( [b -15*b 15*b -b], [-1; 0; 1; 2],  mx, mx );

            params.AI = 1/(12*dx^2)*(Ap-Am);

            row = params.AI(3, 1:5)';

            % upper right block
            params.AI(1, mx-1:mx) = row(1:2);
            params.AI(2, mx)      = row(1);

            % lower left block
            params.AI(mx,   1:2)  = row(4:5);
            params.AI(mx-1, 1)    = row(5);
        end

    elseif( strcmp( params.mbc_type_left,  'dirichlet' ) )

        if( mbc == 4 )

            %%% This section is supposed to provide 6th order accuracy %%%
            Am = spdiags( [-2*b 25*b -245*b 245*b -25*b 2*b], [-3; -2; -1; 0; 1; 2], mx, mx );
            Ap = spdiags( [-2*b 25*b -245*b 245*b -25*b 2*b], [-2; -1;  0; 1; 2; 3], mx, mx );
            params.AI = 1/(180*dx^2)*(Ap-Am);

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



        elseif( mbc == 2 )

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

        else

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Dirichlet and 2nd order Accuracy %%
            Am = spdiags( [ -b b ], [ -1; 0],  mx, mx );
            Ap = spdiags( [ -b b ], [  0; 1],  mx, mx );

            params.AI = 1/(dx^2)*(Ap-Am);

            % upper left block
%           params.AI(1, 1) = (-7/2) / dx^2;
%           params.AI(1, 2) = (-1/2) / dx^2;
%           params.AI(1, 3) = ( 1  ) / dx^2;

%           % lower left block
%           params.AI(mx, mx  ) = (-7/2) / dx^2;
%           params.AI(mx, mx-1) = (-1/2) / dx^2;
%           params.AI(mx, mx-2) = ( 1  ) / dx^2;


            % upper left block
            params.AI(1, 1) = (-9/2) / dx^2;
            params.AI(1, 2) = ( 3/2) / dx^2;

            % lower left block
            params.AI(mx, mx  ) = (-9/2) / dx^2;
            params.AI(mx, mx-1) = ( 3/2) / dx^2;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end

    end

    f = params.eps*params.AI*q;
    f = reshape(f, mx*meqn, 1 );

end
