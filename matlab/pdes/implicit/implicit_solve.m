function y = implicit_solve( yexp, ts, as )
    % Function providing the implicit solve for:
    % y = yexp + dt A_{ii} ( fI( y, t ) );

    global params;
    mx = params.mx;

    y = ( speye(mx) * ( 1 + as/params.eps ) ) \ yexp;

end
