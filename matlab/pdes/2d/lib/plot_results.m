
        % plot the solution
        xc = nodex( (1:mx) + mbc ) + dx/2;
        yc = nodey( (1:my) + mbc ) + dy/2;

size(xc)
size(yc)

        figure(1);
        clf;
%       surf( xc, yc, q );
        figure(2);
        clf;
        contourf( xc, yc, q0 );
        hold off;

%       plot(xc, q(:,1), 'go');
%       hold on;
%       plot(xc, q0, 'r-' );
%       hold on;
%       plot(xc, qex, 'b-' );
%       hold off;
