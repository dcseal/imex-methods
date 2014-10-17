figure(1);
clf;

q1 = q(:,1);
q2 = q(:,2);

plot( tvec, q1, 'go' );
hold on;
plot( tvec, q2, 'ro' );
hold off;

t = ['Van der Pol.  tfinal = ', num2str( tfinal ), '  \epsilon = ', ...
     num2str( params.eps, '%2.1e') ];

title( t );

axis([0 tfinal -5 3]);
