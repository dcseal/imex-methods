function [EvalCoeffs Integrate_Coeffs] = New_Res_Coeffs( spts, local_spts )

    s    = length( spts );
    sloc = length( local_spts );

    EvalCoeffs       = zeros( sloc-1, s );
    Integrate_Coeffs = zeros( sloc-1, s );

    ts = local_spts(1);
    for j=1:(s)

        y = zeros(s,1)';
        y(j) = 1;
        p  = polyfit( spts, y, s-1 );
        pI = polyint( p, 0 );

        ps = polyval( pI, ts );
        for i=1:(sloc-1)

            % todo - can replace one of these EvalCoeffs at time local_spts(end) 
            dt = local_spts(i+1)-local_spts(i);
            Integrate_Coeffs(i,j) = ( ...
                    polyval(pI,local_spts(i+1)) - ps );
            EvalCoeffs(i,j)       = polyval(p, local_spts(i+1));

        end

    end

end
