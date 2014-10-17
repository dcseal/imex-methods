% User supplied Initial condition function
function a = auxinit(x)

    width = 0.2;
    I = find( x > 0.3 & x < 0.5 );
    a = zeros( size(x) );
    a(I) = 1 - sin( pi*( x(I)-0.4 ) / (width) ) .^6 ;

end
