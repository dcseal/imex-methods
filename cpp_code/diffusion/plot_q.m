clear

outputdir = './output/';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial conditions %%
fname = [ outputdir, 'q0000.dat'];

fd = fopen( fname, 'r' );
if( fd == -1 )
    disp(' error opening file');
end
data = fscanf( fd, '%e', [1,Inf] );
fclose(fd);

t0 = data(1);
q0 = data(2:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% final conditions %%
fname = [ outputdir, 'q0001.dat'];

fd = fopen( fname, 'r' );
if( fd == -1 )
    disp(' error opening file');
end

data = fscanf( fd, '%e', [1,Inf] );

fclose(fd);

t = data(1);
q = data(2:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xc = zeros( size(q) );

xlow  = 0.0;
xhigh = 1.0;

sz = size(q);
mx = sz(2);
dx = (xhigh-xlow) / mx;
xc = linspace( xlow, xhigh, mx ) + dx/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the results
figure(1);
clf;
plot(xc, q, 'go' );
hold on;

plot(xc, q0, 'b-');

hold off;

er = norm( q - q0, 1 ) / norm( q0, 1 );
disp( [['Error from initial conditions = ', num2str(er, '%2.5e') ]] );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
