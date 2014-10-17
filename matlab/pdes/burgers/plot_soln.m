% plot the solution
xc = node( (1:mx) + mbc ) - dx/2;

figure(1);
clf;
plot(xc, q(:,1), 'go');
hold on;
plot(xc, qex, 'b-' );
hold on;
plot(xc, q0, 'r-' );
hold off;

title([['Burgers Equation. t = ', num2str(tfinal, '%2.3f'),  ...
'  \epsilon = ', num2str(params.eps, '%2.1e') ]] );
name = [['burgers_eqn_t_', num2str(tfinal, '%2.2f'), '.pdf']];

print(1, name, '-dpdf');
