sp = importdata('point_stress.txt')

subplot(3,1,1)
plot(sp(:,1),sp(:,2))
xlabel('N')
ylabel({'$\sigma_{\varphi}$'},'Interpreter','latex')
title('Circumferential stress $\sigma_{\varphi}$ at point A','Interpreter','latex')

subplot(3,1,2)
plot(sp(:,1),sp(:,3))
xlabel('N')
ylabel({'$\sigma_{\varphi}$'},'Interpreter','latex')
title('Circumferential stress $\sigma_{\varphi}$ at point B','Interpreter','latex')

subplot(3,1,3)
plot(sp(:,1),sp(:,4))
xlabel('N')
ylabel({'$\sigma_{\varphi}$'},'Interpreter','latex')
title('Circumferential stress $\sigma_{\varphi}$ at point C','Interpreter','latex')