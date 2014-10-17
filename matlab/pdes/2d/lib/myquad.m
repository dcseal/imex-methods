function [q, xpts, wgts] = myquad( func, xlow, xhigh, ylow, yhigh )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Gaussian quadrature based on a tensor product of 1D quadrature rules.  The
% parameters are as follows:
%
%    [ xlow(1), xlow(2)   ] = lower rectangular part of cell.
%    [ xhigh(1), xhigh(2) ] = upper right corner of rectangular cell
%
% The output of this function produces the averages value integral:
%
%     1 / dA \iint func(x,y)\ dx dy
%
% Where dA = (xhigh-xlow)*(yhigh-ylow);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [w1d, x1d] = set_gauss_pts();
    [mpts, JUNK] = size(x1d);
    xpts = zeros( mpts^2, 2 );
    wgts = zeros( mpts^2, 1 );

    xc = 0.5*(xlow+xhigh);
    yc = 0.5*(ylow+yhigh);

    % set quadrature weights and points
    l = 1;
    for n=1:mpts
    for k=1:mpts

        xpts(l,1) = xc + 0.5*(xhigh-xlow)*x1d(n);
        xpts(l,2) = yc + 0.5*(yhigh-ylow)*x1d(k);
        wgts(l)   = w1d(n)*w1d(k);
        l = l+1;

    end
    end

    % evaluate user defined funciton at each point
    qpts = func( xpts )';

    % gaussian quadrature.
    q  = 0.25 * qpts * wgts;
     
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
