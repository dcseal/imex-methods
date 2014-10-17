        % plot the solution
        xc = node( (1:mx) + mbc ) - dx/2;

        figure(1);
        clf;

        % numerically computed solution
        plot(xc, q(:,1), 'go');
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

        title([['Scalar Reaction Equation. t = ', num2str(tfinal, '%2.3f') ]] );
        name = [['scalar_reaction_eqn_t_', num2str(tfinal, '%2.2f'), '.pdf']];

%       axis([xlow xhigh -1.1 1.1]);
%       axis([xlow xhigh 0.8 2.1]);
%       print(1, name, '-dpdf');


