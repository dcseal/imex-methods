function y = implicit_solve( yexp, ts, as )
    % Function providing the implicit solve for:
    % y = yexp + dt A_{ii} ( fI( y, t ) );

    global params;
    mx = params.mx;
    my = params.my;

    y = ( speye(mx*my) - as*params.eps*params.AI ) \ yexp;

end
