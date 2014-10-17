% Flux function for the problem
function F = fluxfunc( x, q )

    global params

    F = zeros( size(q) );

%   F(:,1) = params.u * q(:,2);
%   F(:,2) = params.u * q(:,1);

    F(:,1) = params.u * q(:,1);
    F(:,2) = params.u * q(:,2);

end
