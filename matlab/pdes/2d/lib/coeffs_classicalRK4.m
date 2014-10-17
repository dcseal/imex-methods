%%%% Classical RK4 Butcher Tableau %%%%

s = 4;

bE = zeros(s,1);
bE(1) = 1/6;
bE(2) = 1/3;
bE(3) = 1/3;
bE(4) = 1/6;

c = zeros(s,1);
c(2) = 0.5;
c(3) = 0.5;
c(4) = 1;

AE = zeros(s,s);
AE(2,1) = 1/2;

AE(3,2) = 1/2;

AE(4,3) = 1;


% implicit matrices are all zero
bI = zeros(s,1);
AI = zeros( s,s );
