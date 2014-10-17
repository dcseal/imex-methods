function coeffs = Res_Coeffs( spts )
% This function computes the residual coefficients required for integration of the
% polynomial fit to each interval.
%
% Usage:
%
% [coeffs] = Res_Coeffs( spts );
%
% spts = some spacing of points.
%
% coeffs - multiply this by some function values to get an exact integral for
% the polynomial passing through spts.

    s = length( spts );
    coeffs = zeros(s-1,s);
    for i=1:(s-1)

        for j=1:(s)

            y = zeros(s,1)';
            y(j) = 1;

            p = polyfit( spts, y, s-1 );
            pI = polyint( p, 0 );
            coeffs(i,j) = polyval(pI,spts(i+1)) - polyval(pI,spts(i));

        end

    end

end
