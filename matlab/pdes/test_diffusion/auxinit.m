% User supplied Initial condition function
function a = auxinit(x)

%   width = 0.2;
%   I = find( x > 0.3 & x < 0.5 );
%   a = zeros( size(x) );
%   a(I) = 1 - sin( pi*( x(I)-0.4 ) / (width) ) .^6 ;

    global params;
    epsilon = params.eps;

%   a = exp(x).*(1-2*x+3*epsilon*x+epsilon*x.^2);

t = 0;
%a = exp(-t)*(-3*x+x.^2+1+2*epsilon);
a = (-x+x.^2+2*epsilon);

end
