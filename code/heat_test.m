close all; clear all;
figure(1);
heat(1,40,3);

m = 1.6;
s = 0.25;

figure(2);
Q = heat(2,1e2,3);

figure(3);
dim = 16; level = 7;
xi = 0.5*ones(1,dim);
ops = load(['operatorsBlocks',num2str(dim),'_level', num2str(level), '.mat']);  %load ops
ellip2DAffine_heat(xi,ops.ACell,ops.f,level);
xticks(linspace(0,1,5)); yticks(linspace(0,1,5));
title(['Solution u, d=',num2str(length(xi)),', all $\xi_{i} = 0.5$'],'Interpreter','latex');
set(gca,'FontSize',16);

figure(4);
dim = 16; level = 7;
% xi = rand(1,dim)*(10-1) +1;
xi = exp(mvnrnd(m*ones(1,dim), (s^2)*eye(dim)));
ops = load(['operatorsBlocks',num2str(dim),'_level', num2str(level), '.mat']);  %load ops
ellip2DAffine_heat(xi,ops.ACell,ops.f,level);
xticks(linspace(0,1,5)); yticks(linspace(0,1,5));
title(['Solution u, d=',num2str(length(xi)),', $\xi_{i} \sim$ lognormal'],'Interpreter','latex');
set(gca,'FontSize',16);