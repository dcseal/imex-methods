function d = implicit_solve_for_delta( de, q, ts, as )
    % Function providing the implicit solve for:
    % y = yexp + dt A_{ii} ( fI( y, t ) );

    global params;

    d = zeros( size(de) );
    d(2) = de(2) + as;

    q = q(1);
    de = de(1);

    d(1) = ( -(2*q+1/as) + sqrt( (2*q+1/as)^2 + 4*de/as ) ) / (4*q+2/as);

end
