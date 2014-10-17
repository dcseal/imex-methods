        % plot the solution
        xc = node( (1:mx) + mbc ) - dx/2;

        figure(1);
        clf;

        % numerically computed solution
        plot(xc, q, 'go');
        hold on;

        % exact solution (if we have one)
        plot(xc, qex, 'b-' );
        hold on;

        % initial condition (always red)
        plot(xc, q0, 'r-' );
        hold on;

        % plot x-axis
        plot(xc, [0], 'k-' );
        hold off;

        title([['Scalar Advection Equation. t = ', num2str(tfinal, '%2.3f') ]] );

        legend('Approximate Soln', 'Exact Soln', 'Initial Condition');
