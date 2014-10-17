% Explicit right hand side function
function f = fE(t,y)

    global params
    mx   = params.mx;
    my   = params.my;
    meqn = params.meqn;

    % source term gets treated implicitly.
    y = reshape( y, mx, my, meqn );
    f = ConstructL( y, t);
    f = reshape( f, my*mx*meqn, 1 );

end
