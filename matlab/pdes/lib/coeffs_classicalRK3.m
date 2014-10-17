%%%% Classical RK4 Butcher Tableau %%%%

s = 3;

bE = zeros(s,1);
bE(1) = 1/6;
bE(2) = 2/3;
bE(3) = 1/6;

c = zeros(s,1);
c(2) = 0.5;
c(3) = 1.0;

AE = zeros(s,s);
AE(2,1) = 1/2;

AE(3,1) = -1;
AE(3,2) = 2;


% implicit matrices are all zero
bI = zeros(s,1);
AI = zeros( s,s );
