
name = '_imex_eps-1.0e-06';

framenum = 3;
outputdir = ['output' , name ];
meqn = 2;


% this variable is used for read_q()
output_num = 8

% Read in the current value of q
fname = [outputdir , '/q0008.error.dat'];

fd = fopen( fname, 'r' );
mt = fscanf( fd,'%d', [1,1] );
er = fscanf( fd, '%e', [1,inf] );
fclose(fd);

er = reshape( er, [2,mt] )';

tvec = linspace(0, 0.5, mt );
figure(1);
clf;
plot( tvec, er(:,1), 'rx' );
figure(2);
plot( tvec, er(:,2), 'b+' );
