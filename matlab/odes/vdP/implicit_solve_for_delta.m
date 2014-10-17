function d = implicit_solve_for_delta( de, q, ts, as )
    % Function providing the implicit solve for:
    % y = yexp + dt A_{ii} ( fI( y, t ) );

    global params;

    d = de;

    yn = de(1);  z = de(2);
    x  = q(1);  w = q(2);

    rad = ...
    params.eps^2+2*params.eps*as+4*params.eps*as*x*z+4*params.eps*as*x*w+as^2+4*as^2*x*z+4*as^2*x*w+4*as^2*x^2*z*w+4*as^2*x^2*w^2+4*as*z*yn*params.eps+4*as^2*z^2+4*as*w*yn*params.eps+4*as^2*w*z;

    d = (1/2)*(-params.eps-as-2*as*x*z-2*as*x*w+ sqrt( max(rad, 0 ) ))/(as*(z+w));


end
