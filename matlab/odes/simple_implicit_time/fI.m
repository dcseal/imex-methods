% Implicit right hand side function
function f = fI(t,q)

    f = zeros( size(q) );
    f(1) = cos(t);
    f(2) = 1.0;

end
