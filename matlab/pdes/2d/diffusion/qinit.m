% User supplied Initial condition function
function [q] = qinit(xpts)

    [mpts, mdim] = size(xpts);
    q = zeros( mpts, 1 );
    for n=1:mpts
        x = xpts(n,1);
        y = xpts(n,2);
        q(n,1) = sin(2*pi*x)*sin(2*pi*y);
    end


end
