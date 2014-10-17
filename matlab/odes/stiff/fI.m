% Implicit right hand side function
function f = fI(t,q)

    global params;


    f = - q / params.tau;

end
