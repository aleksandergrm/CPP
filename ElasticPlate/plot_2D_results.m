close all;

% plot contour data
N = '50000';
path = strcat('data/N=',N,'/');
X = importdata(strcat(path,'x.dat'));
Y = importdata(strcat(path,'y.dat'));
sr = importdata(strcat(path,'sr.dat'));
sf = importdata(strcat(path,'sf.dat'));
srf = importdata(strcat(path,'srf.dat'));
sx = importdata(strcat(path,'sx.dat'));
sy = importdata(strcat(path,'sy.dat'));
sxy = importdata(strcat(path,'sxy.dat'));

ll = importdata(strcat(path,'results.csv'));

v = linspace(-5,5,10);
db = [4 4; -4 4; -4 -4; 4 -4; 4 4]

figure
plot(ll(:,1),ll(:,2),'.')
hold on
plot(db(:,1),db(:,2),'-b','LineWidth',2)
axis equal

%figure
%contourf(X,Y,sr,linspace(0,1,10))
%colorbar
%title('Radial stress')
%axis equal

%figure
%contourf(X,Y,sf,linspace(-1,1,20))
%colorbar
%title('Circumferential stress')
%axis equal

figure
contourf(X,Y,srf,linspace(-1,1,20))
colorbar
title('Shear stress')
axis equal

figure
contourf(X,Y,sx,20)
hold on
plot(db(:,1),db(:,2),'-b','LineWidth',2)
colorbar
title('X stress')
axis equal

figure
contourf(X,Y,sy,20)
hold on
plot(db(:,1),db(:,2),'-b','LineWidth',2)
colorbar
title('Y stress')
axis equal

figure
contourf(X,Y,sxy,20)
hold on
plot(db(:,1),db(:,2),'-b','LineWidth',2)
colorbar
title('Shear stress - XY')
axis equal