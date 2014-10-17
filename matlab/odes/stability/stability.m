% Copywright: James Rossmanith - University of Wisconsin - Madison

% Stability regions of various methods
clear;

N=100;
xt=linspace(-6.5,1.5,N/2);
yt=linspace(-8,8,N);

[x,y]=meshgrid(xt,yt);

z = x + sqrt(-1)*y;

% Euler
Ax = real(1+z);
Ay = imag(1+z);
Am1 = (Ax.^2 + Ay.^2);

% RK2
Ax = real(1+z+0.5*z.^2);
Ay = imag(1+z+0.5*z.^2);
Am2 = (Ax.^2 + Ay.^2);

% RK3
Ax = real(1+z+0.5*z.^2+z.^3/6);
Ay = imag(1+z+0.5*z.^2+z.^3/6);
Am3 = (Ax.^2 + Ay.^2);

% RK4
Ax = real(1+z+0.5*z.^2+z.^3/6+z.^4/24);
Ay = imag(1+z+0.5*z.^2+z.^3/6+z.^4/24);
Am4 = (Ax.^2 + Ay.^2);

% RK5
Ax = real(1+z+0.5*z.^2+z.^3/6+z.^4/24+z.^5/120);
Ay = imag(1+z+0.5*z.^2+z.^3/6+z.^4/24+z.^5/120);
Am5 = (Ax.^2 + Ay.^2);


%%%%%%%%%%%%%%%%%%%%%
% SDC 2
%%%%%%%%%%%%%%%%%%%%%
Ustar = (1+z);
err = 1 - Ustar + (z./2) .* (1 + Ustar);
Ax = real(Ustar + err);
Ay = imag(Ustar + err);
Asdc2 = (Ax.^2 + Ay.^2);


%%%%%%%%%%%%%%%%%%%%%
% SDC 3
%%%%%%%%%%%%%%%%%%%%%
Morder = 3;

z = (x + sqrt(-1)*y)./2;
u(:,:,1) = ones(size(z));
err(:,:,1) = zeros(size(z));
for n=1:(Morder-1)
    u(:,:,n+1) = (1+z).*u(:,:,n);
end

for k=1:(Morder)
    
    IU(:,:,1) = (z./12).*( -u(:,:,3) + 8.*u(:,:,2) + 5.*u(:,:,1) );
    
    IU(:,:,2) = (z./12).*( -u(:,:,1) + 8.*u(:,:,2) + 5.*u(:,:,3) );

    for n=1:(Morder-1)
        err(:,:,n+1) = (1+z).*err(:,:,n) - (u(:,:,n+1)-u(:,:,n)) ...
            + IU(:,:,n);
    end

    u = u + err;
    
end

Ax = real( u(:,:,Morder) );
Ay = imag( u(:,:,Morder) );
Asdc3 = (Ax.^2 + Ay.^2);


%%%%%%%%%%%%%%%%%%%%%
% SDC 4
%%%%%%%%%%%%%%%%%%%%%
Morder = 4;

z = (x + sqrt(-1)*y)./3;

zz(:,:,1) = 3/4 * z;
zz(:,:,2) = 6/4 * z;
zz(:,:,3) = 3/4 * z;

u(:,:,1) = ones(size(z));
err(:,:,1) = zeros(size(z));
for n=1:(Morder-1)
    u(:,:,n+1) = ( 1 + zz(:,:,n) ).*u(:,:,n);
end

for k=1:(Morder-1)
            
    IU(:,:,1) =  (1/576).*z.*(59*u(:,:,1) + 94*u(:,:,2) + 5*u(:,:,4) - 14*u(:,:,3));
    IU(:,:,2) = -(1/36).*z.*(2*u(:,:,1) - 11*u(:,:,2) + 2*u(:,:,4) - 11*u(:,:,3));
    IU(:,:,3) =  (1/576).*z.*(5*u(:,:,1) - 14*u(:,:,2) + 59*u(:,:,4) + 94*u(:,:,3));

    for n=1:(Morder-1)
        err(:,:,n+1) = ( 1 + zz(:,:,n) ).*err(:,:,n) - (u(:,:,n+1)-u(:,:,n)) ...
            + IU(:,:,n);
    end

    u = u + err;
    
end

Ax = real( u(:,:,Morder) );
Ay = imag( u(:,:,Morder) );
Asdc4 = (Ax.^2 + Ay.^2);


%%%%%%%%%%%%%%%%%%%%%
% SDC 5
%%%%%%%%%%%%%%%%%%%%%
Morder = 5;

z = (x + sqrt(-1)*y)./4;

zz(:,:,1) = (2-sqrt(2))*z;
zz(:,:,2) = sqrt(2)*z;
zz(:,:,3) = sqrt(2)*z;
zz(:,:,4) = (2-sqrt(2))*z;

u(:,:,1) = ones(size(z));
err(:,:,1) = zeros(size(z));
for n=1:(Morder-1)
    u(:,:,n+1) = ( 1 + zz(:,:,n) ).*u(:,:,n);
end

for k=1:(Morder-1)                              
    
    IU(:,:,1) = (z./480).*((23+4*sqrt(2))*u(:,:,1)+(-13*sqrt(2)+64)*u(:,:,2)+...
        (-72*sqrt(2)+96)*u(:,:,3)+(-43*sqrt(2)+64)*u(:,:,4)+...
        (-7+4*sqrt(2))*u(:,:,5));

    IU(:,:,2) = (sqrt(2)*z./960).*((-8-15*sqrt(2))*u(:,:,1)+146*u(:,:,2)+...
        144*u(:,:,3)-34*u(:,:,4)+...
        (-8+15*sqrt(2))*u(:,:,5));

    IU(:,:,3) = (sqrt(2)*z./960).*((-8+15*sqrt(2))*u(:,:,1)-34*u(:,:,2)+...
        144*u(:,:,3)+146*u(:,:,4)+...
        (-8-15*sqrt(2))*u(:,:,5));

    IU(:,:,4) = (z./480).*((-7+4*sqrt(2))*u(:,:,1)+(-43*sqrt(2)+64)*u(:,:,2)+...
        (-72*sqrt(2)+96)*u(:,:,3)+(-13*sqrt(2)+64)*u(:,:,4)+...
        (23+4*sqrt(2))*u(:,:,5));
                              

    for n=1:(Morder-1)
        err(:,:,n+1) = ( 1 + zz(:,:,n) ).*err(:,:,n) - (u(:,:,n+1)-u(:,:,n)) ...
            + IU(:,:,n);
    end

    u = u + err;
    
end

Ax = real( u(:,:,Morder) );
Ay = imag( u(:,:,Morder) );
Asdc5 = (Ax.^2 + Ay.^2);


%figure(1)
%clf
%contour(x,y,Am1,[1 1],'r');
%axis('equal');
%axis([-6.5 1.5 -8 8]);
%hold on;
%contour(x,y,Am2,[1 1],'g');
%contour(x,y,Am3,[1 1],'b');
%contour(x,y,Am4,[1 1],'r');
%contour(x,y,Am5,[1 1],'g');
%hold off;
%axis on; grid on; box on;
%set(gca,'fontsize',16);
%set(gca,'xtick',-9:1:3);
%set(gca,'ytick',-8:1:8);

%figure(2)
%clf
%contour(x,y,Am1,[1 1],'r');
%axis('equal');
%axis([-6.5 1.5 -8 8]);
%hold on;
%contour(x,y,Asdc2,[1 1],'g');
%contour(x,y,Asdc3,[1 1],'b');
%contour(x,y,Asdc4,[1 1],'r');
%contour(x,y,Asdc5,[1 1],'g');
%hold off;
%axis on; grid on; box on;
%set(gca,'fontsize',16);
%set(gca,'xtick',-9:1:3);
%set(gca,'ytick',-8:1:8);

%figure(3)
%clf
for k=1:4
    %figure(3)
    %subplot(1,4,k);
    figure(k)
    plot([-6.5 1.5],[0 0],'k--');
    hold on;
    plot([0 0],[-8 8],'k--');
    if k==1
        [tmp,c1]=contour(x,y,Am2,[1 1],'r--');
        hold on;
        [tmp,c2]=contour(x,y,Asdc2,[1 1],'b');
        t1=title('Second Order Method');
    elseif k==2
        [tmp,c1]=contour(x,y,Am3,[1 1],'r--');
        hold on;
        [tmp,c2]=contour(x,y,Asdc3,[1 1],'b');
        t1=title('Third Order Method');
    elseif k==3
        [tmp,c1]=contour(x,y,Am4,[1 1],'r--');
        hold on;
        [tmp,c2]=contour(x,y,Asdc4,[1 1],'b');
        t1=title('Fourth Order Method');
    else
        [tmp,c1]=contour(x,y,Am5,[1 1],'r--');
        hold on;
        [tmp,c2]=contour(x,y,Asdc5,[1 1],'b');
        t1=title('Fifth Order Method');
    end   
    set(c1,'linewidth',2);
    set(c2,'linewidth',2);
    set(t1,'fontsize',16);
    hold off;
    axis('equal');
    axis([-6.5 1.5 -8 8]);
    axis on; grid off; box on;
    set(gca,'fontsize',16); 
    set(gca,'xtick',-9:1:3);    
    set(gca,'ytick',-8:1:8);
end

figure(1);
