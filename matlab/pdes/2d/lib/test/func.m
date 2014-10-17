function q = func(xpts)

    [mpts, md] = size(xpts);
    q = zeros( mpts, 1 );
    for n=1:mpts
        x = xpts(n,1);
        y = xpts(n,2);
        q(n,1) = y*sin(x) + x * cos(y);
        q(n,2) = cos(x);
    end

end
