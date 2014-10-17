mx = 10;
my = 6;

G = ones( mx+1, my+1, 2);
F = ones( size(G) );

I = 1:mx;  J = 1:my;

G(I+1,J,1) - G(I,J,1)
