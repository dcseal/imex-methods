
    % plot the solution
    xc = node( (1:mx) + mbc ) - dx/2;

    figure(1);
    clf;
    plot(xc, q0, 'r-' , 'Linewidth', 3);
    hold on;
    plot(xc, q(:,1), 'go', 'Linewidth', 3 );
%   hold on;
%   plot(xc, qex, 'b-' );
    hold off;

    axis([xlow xhigh -0.1 1.1]);

    l = legend('Initial Condition', 'Final Condition');
    
    t = title([['Convection-Diffusion Equation.  t = ', num2str(tfinal,'%2.2f')]] );
