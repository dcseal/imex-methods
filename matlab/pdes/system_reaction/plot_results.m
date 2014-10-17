        % plot the solution
        xc = node( (1:mx) + mbc ) - dx/2;

        figure(1);
        clf;

        % numerically computed solution
        plot(xc, q(:,1), 'bo', 'Linewidth', 2);
        hold on;
        plot(xc, q(:,2), 'b+', 'Linewidth', 2);

        % initial condition (always red)
        plot(xc, q0, 'r-' , 'Linewidth', 3);
        hold on;


        % exact solution (if we have one)
%       plot(xc, qex, 'b-' );
%       hold on;

        % plot x-axis
        plot(xc, [0], 'k-' );
        hold off;

        title([['Systems Reaction Equation. t = ', num2str(tfinal, '%2.3f') ]] );
        name = [['system_reaction_eqn_t_', num2str(tfinal, '%2.2f'), '.eps']];

        legend('u', 'v', 'Initial Condition' );
        %legend('u', 'v', 'Exact Soln', '', 'Initial Condition');

        axis([xlow xhigh -1.1 1.1]);
        print(1, '-depsc', name);


