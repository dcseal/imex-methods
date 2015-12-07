function [spts s] = get_quad_pts(N)
% This function determines the quadrature points used for integration of the
% residual.  The default behaviour is to use uniform time points.

    if( N == 4 )

        % These weights and points integrate a polynomial of degree 5 exactly,
        % and therefore are 6th order accurate.
%       spts = [ -1.00000000000000000000e+00, -4.47213595499957927704e-01, ...
%                4.47213595499957927704e-01, 1.00000000000000000000e+00];
%       spts = [-1, -1/sqrt(3), 1/sqrt(3), 1];
%       spts = [-1 -0.5 0.5 1];
        spts = linspace(-1,1,4);

        wgts = [1.66666666666666657415e-01, 8.33333333333333370341e-01, ...
                8.33333333333333370341e-01, 1.66666666666666657415e-01];

    elseif( N == 3 )

        % These weights and points integrate a polynomial of degree 3 exactly,
        % and therefore are 4th order accurate.
        spts = [-1.0, 0.0, 1.0];

        wgts = [ 3.33333333333333314830e-01, ...
                 1.33333333333333325932e+00, ...
                 3.33333333333333314830e-01];

    elseif( N == 2 )
    
        % 2nd Order SDC
        spts = [-1 1];
        wgts = [1 1];

    else
        spts = linspace(-1,1,N);
    end

    s = length( spts );

end
