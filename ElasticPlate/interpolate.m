close all;

X = importdata('x.dat');
Y = importdata('y.dat');
SX = importdata('sx.dat');
SY = importdata('sy.dat');
SXY = importdata('sxy.dat');

[xq,yq] = meshgrid(-2:.1:2, -2:.1:2);
SXq = griddata(X,Y,SX,xq,yq,'v4');

figure
mesh(xq,yq,SXq);