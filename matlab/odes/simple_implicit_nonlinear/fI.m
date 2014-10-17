% Implicit right hand side function
function f = fI(t,q)

    f = zeros( size(q) );
    f(1) = -q(1)^2;
    f(2) = 1.0;

end
