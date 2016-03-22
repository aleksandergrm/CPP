figure
M1 = importdata('data/line_data/line_results_03.dat');
M2 = importdata('data/line_data/line_results_07.dat');
M3 = importdata('data/line_data/line_results_10.dat');
M4 = importdata('data/line_data/line_results_50.dat');
M5 = importdata('data/line_data/line_results_100.dat');
M6 = importdata('data/line_data/line_results_500.dat');
M7 = importdata('data/line_data/line_results_50000.dat');

plot(M1(:,1),M1(:,4),M2(:,1),M2(:,4),M3(:,1),M3(:,4),M4(:,1),M4(:,4),M5(:,1),M5(:,4),M6(:,1),M6(:,4),M7(:,1),M7(:,4))
axis([1.0 1.1 0.1 -inf])
xlabel({'$\rho$'},'Interpreter','latex')
ylabel({'$\tau_{\rho\varphi}$'},'Interpreter','latex')
legend('N=3','N=7','N=10','N=50','N=100','N=500','N=50000')
title('Shear stress - $\tau_{\rho\varphi}$ (point B: $\rho\in[1.0,1.2]$, $\varphi=\frac{\pi}{4}$','Interpreter','latex')