% initial condition:


% method of Hochbruck et a. (Exp4)
% M. Hochbruck, C. Lubich, H. Selhofer, ``Exponental Integrators for large
% systems of differential equations'', SIAM Journal of Scientific Computing, 19,
% (1998)

% phi1(z) = ( exp(z) - 1 ) / z = 1 + z / 2! + z^2 / 3! + z^3 / 4! + z^4 / 5! + \cdots

f0 = frhs( y0 );

At = 1/3*h*A0;
phi1 = 1 + At / 2 + At^2 / 6 + At^3 / 24 + At^4 / 120;
k1 = At * f0;

At = 2/3*h*A0;
phi2 = 1 + At / 2 + At^2 / 6 + At^3 / 24 + At^4 / 120;
k2 = phi;

At = h*A0;
phi3 = 1 + At / 2 + At^2 / 6 + At^3 / 24 + At^4 / 120;
k2 = phi;

w4 = -7/300*k1 + 97/150*k2 - 37/300*k3;
u4 = y0 + h*w4;
r4 = frhs(u4) - f0 - h*A0*w4;

k4 = phi1*r4; k5 = phi2*r4; k6 = phi3*r4;

w7 = 59/300*k1 - 7/75*k2 + 269/300*k3 + 2/3*( k4+k5+k6 );
u7 = y0 + h*w7;
r7 = frhs(y7) - f0 - h*A0*w7;

k7 = phi1*r7;

y = y0 + h*(k3+k4-4/3*k5 + k6 + k7/6 );
