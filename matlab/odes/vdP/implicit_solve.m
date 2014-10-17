function y = implicit_solve( yexp, ts, as )
    % Function providing the implicit solve for:
    % y = yexp + dt A_{ii} ( fI( y, t ) );

    global params;

    y = yexp;
    y(2) = ( params.eps * yexp(2) - as*yexp(1) ) / ...
    ( params.eps         - as*(1-yexp(1)^2) );

end
