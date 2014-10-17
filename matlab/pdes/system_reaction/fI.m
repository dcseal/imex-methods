% Implicit right hand side function
function f = fI(t,q)


    global params
    mx   = params.mx;
    aux  = params.aux;
    meqn = params.meqn;

    q = reshape(q, mx, meqn );

    f = zeros( size(q) );
%   f(:,1) = zeros( size(q(:,1) ) );
%   f(:,2) = -(q(:,2) - params.r * q(:,1) ) / params.tau;
    f = reshape( f, mx*meqn, 1 );


end
