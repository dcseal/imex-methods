clear all;

global params

params.mbc    = 1;
params.meqn   = 1;
params.sorder = 2;
params.mx     = 2;
params.my     = 3;


addpath('../');
A = [1 2 3; 4 5 6];

Abig = fill_boundary( A, 0 );

disp(['A = ']);
A

disp(['Abig = ']);
Abig
