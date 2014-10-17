mx = 100;

e = ones(mx,1);
A = spdiags([e -2*e e], -1:1, mx, mx );

f = A * e;
u = zeros( mx,1 );



r = f - A*u;
p = r;

tol = 1e-10;
numiter = 0;
while( norm(r,1) > tol )

    w = A*p;
    alpha = (r'*r) / (p'*w);
    u = u + alpha * p;

    rkp1 = r - alpha * w;

    % check error
    beta = (rkp1'*rkp1) / (r'*r);
    p    = rkp1 + beta * p;

    r = rkp1;

    numiter = numiter + 1;

end

err = norm( u - e, 1 );

disp([ ['numiter = ', num2str(numiter,'%d'), '; err = ', num2str(err, '%2.3e')]]);
