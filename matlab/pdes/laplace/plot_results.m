
        % plot the solution
        xc = node( (1:mx) + mbc ) - dx/2;

        figure(1);
        clf;
        plot(xc, q(:,1), 'go');
        hold on;
        plot(xc, q0, 'r-' );
        hold on;
        plot(xc, qex, 'b-' );
        hold off;
