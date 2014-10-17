function q = implicit_solve( qn, ts, as )
    % Function providing the implicit solve for:
    % y = yexp + dt A_{ii} ( fI( y, t ) );

    global params;
    mx   = params.mx;
    tau  = params.tau;
    aux  = params.aux;
    meqn = params.meqn;


    qn = reshape(qn, mx, meqn );

    u = qn(:,1);
    v = qn(:,2);

    q(:,1) = u;
    q(:,2) = ( tau * v + as*params.r * u ) / ( tau+as ) ;

    q = reshape( q, mx*meqn, 1 );

end
