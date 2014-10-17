function y = implicit_solve( yexp, ts, as )

    % implicit solve for:
    % y = yexp + dt A_{ii} ( f(y,ts) )

    y = zeros( size(yexp) );
    y(1) = yexp(1) + as * cos(ts);
    y(2) = yexp(2) + as;

end
