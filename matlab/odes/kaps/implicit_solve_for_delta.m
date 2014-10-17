function d = implicit_solve_for_delta( de, q, ts, as )
    % Function providing the implicit solve for:
    % y = yexp + dt A_{ii} ( fI( y, t ) );

    global params;

    d = de;

    yn = de(1);  z = de(2);
    x  = q(1);  w = q(2);

%   d(1) = ( de(1) + as/params.eps*( (q(2)+de(2))^2 - q(2)^2 ) ) / (1 + as/params.eps );

    d(1) = (yn*params.eps+as*z^2+2*as*z*w)/(params.eps+as);

end
