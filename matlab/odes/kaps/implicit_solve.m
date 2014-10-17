function y = implicit_solve( yexp, ts, as )
    % Function providing the implicit solve for:
    % y = yexp + dt A_{ii} ( fI( y, t ) );

    global params;

    y = yexp;

    y(1) = ( params.eps * yexp(1) + as*yexp(2)^2 ) / ( params.eps + as );

end
