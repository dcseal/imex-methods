% Implicit right hand side function
function f = fI_delta(t,de,q)

    global params

    f = zeros( size(q) );
    f(1) = 1/params.eps*( ( (q(2)+de(2))^2 - q(2)^2 ) - de(1) );
    f(2) = 0;

end
