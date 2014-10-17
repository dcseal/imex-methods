function [IM EM] = polynom_matrices( sorder, num_stages, spts, long_spts )
% This function computes the integration and evaluation matrices necessary
% for evaluating the Lagrange interpolant at each of the sub time points

    IM = zeros( num_stages*(s-1), s);
    EM = zeros( num_stages*s,     s);

    for j=1:(s)

        y = zeros(s,1)';
        y(j) = 1;

        p = polyfit( spts, y, s-1 );
        pI = polyint( p, 0 );

        for i=1:(s-1)

            IM(i,j) = polyval(pI,long_spts(i+1)) - polyval(pI,long_spts(i));
            EM(i,j) = polyval(pI,long_spts(i));

        end

        EM(i,j) = polyval(pI,spts(i));

    end

end
