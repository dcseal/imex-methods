function y = implicit_solve( yexp, ts, as )
    % Function providing the implicit solve for:
    % y = ys + dt * A_{ii} ( fI( y, t ) );

    global params;

    y = yexp;

    y = ( yexp ) / ( 1 + as/params.tau );

end
