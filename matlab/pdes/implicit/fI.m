% Implicit right hand side function
function f = fI(t,q)

    global params
    mx   = params.mx;
    meqn = params.meqn;
    dx   = params.dx;
    mbc  = params.mbc;

    f = -q / params.eps;
    f = reshape(f, mx*meqn, 1 );

end
