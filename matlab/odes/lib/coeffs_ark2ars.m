%%%% ARK2ARS IMEX RK method %%%%%
% Ascher, Ruuth, Spiteri.  `Implicit-explicit Runge-Kutta methods for
% time-dependent partial differential equation.'

s     = 3;

gamma = 1 - 0.5*sqrt(2);
d     = -2*sqrt(2)/3;

bE    = zeros(s,1);
bE(1) = 0;
bE(2) = 1-gamma;
bE(3) = gamma;

% Intermediate time values
c = zeros(size(bE));
c(2) = gamma;
c(3) = 1;


%%% Explicit Matrix %%%
AE = zeros(s,s);

AE(2,1) =  gamma;

AE(3,1) = d;
AE(3,2) = 1-d;

%%% Implicit Matrix %%%
AI = zeros(s,s);
AI(2,2) = gamma;

AI(3,2) = 1-gamma;
AI(3,3) = gamma;

bI = zeros(s,1);
bI(1) = 0;
bI(2) = 1-gamma;
bI(3) = gamma;
