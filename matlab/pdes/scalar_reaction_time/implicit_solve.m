function y = implicit_solve( yexp, ts, as )
    % Function providing the implicit solve for:
    % y = yexp + dt A_{ii} ( fI( y, t ) );

    global params;
    mx  = params.mx;
    aux = params.aux;
    tau = params.tau;

    % First problem (meqn = 1):
    y = ( tau*yexp ) / (tau+as);

end
