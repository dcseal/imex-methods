% Source term function
function psi = Psi(t, q, params )

    u        = params.u;
    dx       = params.dx;
    meqn     = params.meqn;
    mx       = params.mx;
    mbc      = params.mbc;

%   psi      = zeros( size(q) );
%   psi(1,1) = u/dx*sin(2*pi*t);

    % Stiff source term here
    aux = params.aux;
    psi = (-1/params.tau)*( q - aux );

end
