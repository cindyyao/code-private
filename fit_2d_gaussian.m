function params = fit_2d_gaussian(data)
%% ---Fitting Functions---
%
% Coeficients A convention:
%	A = [Amplitude, x0, x-Width, y0, y-Width, Angle(in Radians), b0]
%
% X-data convention:
%	X is of size(n,n,2) where 
%	X(:,:,1) : x-coordinates,  
%	X(:,:,2) : y-coordinates.
%
% 2D Rotated Gaussian function ( A requires 6 coefs ).
f = @(A,X) A(1)*exp( -(...
    ( X(:,:,1)*cos(A(6))-X(:,:,2)*sin(A(6)) - A(2)*cos(A(6))+A(4)*sin(A(6)) ).^2/(2*A(3)^2) + ... 
    ( X(:,:,1)*sin(A(6))+X(:,:,2)*cos(A(6)) - A(2)*sin(A(6))-A(4)*cos(A(6)) ).^2/(2*A(5)^2) ) ) + A(7);

%% Parameters
[n, m] = size(data);
A0(1) = max(data(:)); % initial condition
[~, A0(2)] = max(max(data));
[~, A0(4)] = max(max(data, [], 2));
A0(3) = min(m,n)/4;
A0(5) = min(m,n)/4;
A0(6) = 0;
A0(7) = 0;

%% ---Build numerical Grids---
% Numerical Grid
[x,y]=meshgrid(1:n,1:m); X=zeros(m,n,2); X(:,:,1)=x; X(:,:,2)=y;
% High Resolution Grid
h=3; [xh,yh]=meshgrid(1/h:1/h:n,1/h:1/h:m); Xh=zeros(h*m,h*n,2); Xh(:,:,1)=xh; Xh(:,:,2)=yh;

%% ---Fit---
% Define lower and upper bounds [Amp,xo,wx,yo,wy,fi]
lb = [0,1,0,1,0,0,0];
ub = [realmax('double'),n,n,m,m,pi,realmax('double')];

[A,resnorm,res,flag,output] = lsqcurvefit(f,A0,X,data,lb,ub);

% disp(output); % display summary of LSQ algorithm
params.A = A(1); params.x0 = A(2); params.xd = A(3);
params.y0 = A(4); params.yd = A(5); params.angle = A(6);
params.b0 = A(7);

%% ---Plot Data---
% Plot 3D Data and Fitted curve
% hf1=figure(1); set(hf1,'Position',[1000 600 800 500]); 
% C=del2(f(A,Xh)); mesh(xh,yh,f(A,Xh),C); hold on
% 
% surface(x,y,data,'EdgeColor','none'); alpha(0.5); 
% colormap('pink'); view(-60,20); grid on; hold off


end