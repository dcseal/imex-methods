% Flux function for the problem
function F = fluxfunc( x, q )

    global params

    F = zeros( size(q) );

    F(:,1) = params.a * q(:,2);
    F(:,2) = params.b * q(:,1);

end
