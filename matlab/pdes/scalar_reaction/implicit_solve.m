function y = implicit_solve( yexp, ts, as )
    % Function providing the implicit solve for:
    % y = yexp + dt A_{ii} ( fI( y, t ) );

    global params;
    mx = params.mx;
    tau = params.tau;
    aux = params.aux;

    % y = ( tau*yexp + as*aux ) / (tau+as);
    y = ( tau*yexp ) / (tau+as);

end
