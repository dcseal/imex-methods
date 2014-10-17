function y = implicit_solve( yexp, ts, as )
    % Function providing the implicit solve for:
    % y = yexp + dt A_{ii} ( fI( y, t ) );

    global params;
    mx = params.mx;

    y = ( eye( mx ) - params.uimp * as*params.A ) \  yexp;

end
