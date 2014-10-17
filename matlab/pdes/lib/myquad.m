function [q, w1d, x1d] = myquad( func, a, b, tol )

    [w1d, x1d] = set_gauss_pts();

    xc = 0.5*(a+b);    
    xpts = xc + 0.5*(b-a)*x1d;
    qpts = func( xpts )';
   
    [meqn, mpts] = size(qpts);
    
    q = zeros( meqn, 1 );
    q = 0.5 * (b-a) * qpts * w1d;
     
end


function [w1d, x1d] = set_gauss_pts( mpoints )

    mpoints = 4;
    x1d = zeros(mpoints,1);
    w1d = zeros(mpoints,1);

    sq2 = sqrt(2.0);
    sq3 = sqrt(3);
    sq5 = sqrt(5);
    sq7 = sqrt(7);

    x1d(1) = sqrt( ( 3.0 - 2.0*sq2*sq3/sq5 ) / 7.0 );
    x1d(2) = - x1d(1);

    x1d(3) = sqrt( ( 3.0 + 2.0*sq2*sq3/sq5 ) / 7.0 );
    x1d(4) = - x1d(3);

    w1d(1) = ( 18.0+sqrt( 30.0 ) ) / 36.0;
    w1d(2) = ( 18.0+sqrt( 30.0 ) ) / 36.0;

    w1d(3) = ( 18.0-sqrt( 30.0 ) ) / 36.0;
    w1d(4) = ( 18.0-sqrt( 30.0 ) ) / 36.0;


end
