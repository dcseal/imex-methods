%%%% Classical RK4 Butcher Tableau %%%%

s = 2;

bE = zeros(s,1);
bE(1) = 1/2;
bE(2) = 1/2;

c = zeros(s,1);
c(2) = 1.0;

AE = zeros(s,s);
AE(2,1) = 1.0;


% implicit matrices are all zero
bI = zeros(s,1);
AI = zeros( s,s );
