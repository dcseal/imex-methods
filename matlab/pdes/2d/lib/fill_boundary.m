function [qbig] = fill_boundary( q, t )

    global params

    mbc    = params.mbc;
    meqn   = params.meqn;
    sorder = params.sorder;
    mx     = params.mx;
    my     = params.my;

    qbig   = zeros( mx + 2*mbc, my + 2*mbc, meqn );

    % center elements are identical
    I = 1:mx;
    J = 1:my;
    qbig( I+mbc, J+mbc, : ) = q(I,J,:);

    %%%%%%%% Fill in the ghost cells %%%%%%%%%%%%%%%%
    %if( strcmp( params.mbc_type_left , 'periodic' ) )

        % left and right values
        for j=1:my
        for n=1:mbc

            qbig(n,        mbc+j, :)  = q(mx-mbc+n, j,  : );
            qbig(mx+mbc+n, mbc+j, :)  = q(n,         j,  : );

        end
        end

        % top and bottom values
        for i=1:mx
        for j=1:mbc

            qbig(i+mbc,        j, :)  = q(i, my-mbc+j,  : );
            qbig(i+mbc, my+mbc+j, :)  = q(i,        j,  : );

        end
        end

        % TODO - DEAL WITH SETTING THE CORNER'S

        return;

    %end

end
