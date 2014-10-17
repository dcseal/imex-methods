nout = 5;

xlow   = 0;
ylow   = 0;
xhigh  = 2*pi;
yhigh  = 2*pi;

mxstart = 2;
mystart = 2;

for n=1:nout

    er = zeros(nout,1);

    mx = mxstart*2^(n-1);
    my = mystart*2^(n-1);

    dx = (xhigh-xlow)/mx;
    dy = (yhigh-ylow)/my;

    q   = zeros(mx,my,2);
    for i=1:mx
    for j=1:my

        xl = xlow + dx * (i-1);
        yl = ylow + dy * (j-1);

        q(i,j,:)    = myquad( @func, xl, xl+dx, yl, yl+dx );
        %qex(i,j)  = dblquad(@(x,y)y*sin(x) + x * cos(y), xl, xl+dx, yl, yl+dy ) / (dx*dy);


    end
    end

%   er(n) = norm( qex-q, 2 ) / norm( qex, 2 );
%   if( n > 1 )
%       disp([[' log(ratio) = ', num2str( log2( er(n-1) / er(n) ), '%2.5e' ) ]] );
%   end

end

xc = (1:mx) * dx;
yc = (1:my) * dy;

surf( xc, yc, q(:,:,1) );
surf( xc, yc, q(:,:,2) );
