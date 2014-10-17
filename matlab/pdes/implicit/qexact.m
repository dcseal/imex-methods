% User supplied Initial condition function
function q = qexact(a,b)

    global params
    mx     = params.mx;
    meqn   = params.meqn;

    q = quad( @q2func, a, b );

end

function q1 = q1func(x)
% q = sin( \lambda x )

    global params

    lambda = 2*pi;
    tf = params.tfinal;
    q1 = ( sin( lambda*(x- params.u*tf) ) + ...
            0.01*cos( lambda*(x-params.u*tf) ) ) ...
            .* exp( - params.eps * lambda^2 * tf );

end

function q2 = q2func(x)


    global params
    q2 = exp( -params.tfinal / params.eps ) * ones( size(x) );

end
