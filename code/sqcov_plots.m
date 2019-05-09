close all; clearvars;
addpath('../ellip/');
addpath('../sim data/');

d = 16; %dimension
N = 5e1; %number of samples
t = 0.315; %threshold (see above)
tname = num2str(t); tname(1:2) = [];

load(['sqcovdata_dim',num2str(d),'threshold0point',num2str(tname),'.mat']);
MFCE_sqcovs = sqcov_data.MFCE_sqcovs;
MFCE_avgtimes = sqcov_data.MFCE_avgtimes;
CE_sqcovs = sqcov_data.CE_sqcovs;
CE_avgtimes = sqcov_data.CE_avgtimes;


figure(1);
loglog(CE_avgtimes,CE_sqcovs,'-s','LineWidth',1.5); hold on;
loglog(MFCE_avgtimes,MFCE_sqcovs,'-d','LineWidth',1.5);
% semilogy(CE_avgtimes,CE_sqcovs,'-s','LineWidth',1.5); hold on;
% semilogy(MFCE_avgtimes,MFCE_sqcovs,'-d','LineWidth',1.5);
xlabel('total runtime'); ylabel('SQCoV')
ylim([1e-5, 1e1]); yticks(10.^(-5:1:1))
% xticks([5e0, 1e1, 5e1, 1e2, 5e2]);
legend('high-fidelity CE', 'MFCE');
title(['Coeff of Var.,d=',num2str(d),',t=',num2str(t)])
set(gca,'FontSize',16);
grid on;