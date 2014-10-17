%{
names(1) =  '_imex_eps-1.0e+00';
names(2) =  '_ark32_eps-1.0e+00';
names(3) =  '_ark43_eps-1.0e+00';
names(4) =  '_imex_eps-1.0e-06';
names(5) =  '_ark32_eps-1.0e-06';
names(6) =  '_ark43_eps-1.0e-06';
%}

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

%   n = 0
%   error = np.zeros( (mt, meqn) )
%   for line in fd:
%       line = line.split()
%       for m in range(meqn):
%           error[n, m] = float( line[m] )
%       n = n+1
%   fd.flush()
%   fd.close()

%   tvec = np.linspace(0, 0.5, mt )

%   plt.plot( tvec, error[:,0], 'go' )
%   print error
%   print tvec
%   plt.show()
