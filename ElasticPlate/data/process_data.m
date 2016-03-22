N = [3 7 50 100 500 1000 5000 10000 50000 100000 200000];

fileName = strcat('N=',num2str(N(2)),'/line_results.dat');
M = dlmread(fileName);

close all

figure
plot(M(:,1),M(:,4))

%S = dlmread('stress_at_points.txt');

%figure
%plot(S(:,1),S(:,2))
%figure
%plot(S(:,1),S(:,3))
%figure
%plot(S(:,1),S(:,4))
