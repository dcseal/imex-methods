function y = implicit_solve( yexp, ts, as )

    % implicit solve for:
    % y = yexp + dt A_{ii} ( f(y,ts) )
    y = zeros( size(yexp) );

    y(1) = ( -1+sqrt(1+4*as*yexp(1)) ) / (2*as);
    y(2) = yexp(2) + as;

end
