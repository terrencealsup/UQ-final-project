close all; clearvars;
addpath('../ellip/');
addpath('../sim data/');
%test coeff of variation


d = 1; %dimension
N = 5e1; %number of samples
t = 0.90; %threshold (see above)
tname = num2str(t); tname(1:2) = [];
%change truth based on above
num_times = 30;
load(['CEdim',num2str(d),'threshold0point',tname,'.mat']);

% [c,avg_time] = coeffofvarMFCE(CE_truth,d,N,t,num_times);
% disp(['Squared coefficient of variation estimate is ', num2str(c)]);
% disp(['runtime is ',num2str(avg_time)]);

disp(d)
% d = 1
MFCE_N_vec = [5e1,5e2,5e3,5e4];
CE_N_vec = [1e2,1e3,1e4,1e5];

% % d = 16
% MFCE_N_vec = [5e3,1e4,5e4,1e5];
% CE_N_vec = [5e3,1e4,5e4,1e5];


MFCE_sqcovs = zeros(1,length(MFCE_N_vec)); MFCE_avgtimes = zeros(1,length(MFCE_N_vec));
CE_sqcovs = zeros(1,length(MFCE_N_vec)); CE_avgtimes = zeros(1,length(MFCE_N_vec));

tic;
for i=1:length(MFCE_N_vec)
    [MFCE_sqcovs(i),MFCE_avgtimes(i)] = coeffofvarMFCE(CE_truth,d,MFCE_N_vec(i),t,num_times);
    [CE_sqcovs(i),CE_avgtimes(i)] = coeffofvarCE(CE_truth,d,CE_N_vec(i),t,num_times);
end
timed = toc;

sqcov_data.MFCE_sqcovs = MFCE_sqcovs;
sqcov_data.MFCE_avgtimes = MFCE_avgtimes;
sqcov_data.CE_sqcovs = CE_sqcovs;
sqcov_data.CE_avgtimes = CE_avgtimes;
sqcov_data.d = d;
sqcov_data.t = t;


disp(['computing squared coeff. of var. took t=',num2str(timed),' seconds, using d=',num2str(d),...
    ',t=',num2str(t)]);
% str = input('Do you want to save current data? y / n: ', 's');
str = 'y'; %autosave
if str == 'y'
    save(['sqcovdata_dim',num2str(d),'threshold0point',num2str(tname),'.mat'],'sqcov_data');
    disp('Saved data');
else
    disp('Did not save data');
end


figure(1);
loglog(CE_avgtimes,CE_sqcovs,'-s','LineWidth',1.5); hold on;
loglog(MFCE_avgtimes,MFCE_sqcovs,'-d','LineWidth',1.5);
xlabel('total runtime'); ylabel('squared coeff. of var.')
legend('high-fidelity CE', 'MFCE');
title(['Coeff of Var.,d=',num2str(d),',t=',num2str(t)])
set(gca,'FontSize',16);
grid on;



function vals = f(x,l,d)
    ops = load(['operatorsBlocks',num2str(d),'_level',num2str(l),'.mat']);
    vals = ellip2DAffine(x,ops.ACell,ops.f,l);
end

function [coeffofvar, avg_time] = coeffofvarCE(truth,d,N,t,num_times)
m = 1.6*ones(1,d); s = 0.25;
if d > 1
    s = (s^2)*eye(d);
end
if d == 1
    p = @(x) normpdf(log(x), m, s) ./ prod(x,2); 
else
    p = @(x) mvnpdf(log(x), m, s) ./ prod(x,2);
end

L = 5; % highest level
rho = 0.9; % technical term 1
dt = 0.01; % technical term 2
ellipfn = @(x) f(x,L,d); % model operator at highest level, current dim

estimates = zeros(1,num_times);
runtimes = zeros(1,num_times);
for i=1:num_times
   Z = exp(mvnrnd(m,s,N)); %initial log-normal sample 
   tic; 
   [estimates(i),~,~] = CE_LogNormal(Z, t, p, ellipfn, N, rho, dt,{});
   runtimes(i) = toc;
end
coeffofvar = mean((estimates - truth).^2) / (truth^2);
avg_time = mean(runtimes);
end

function [coeffofvar, avg_time] = coeffofvarMFCE(truth,d,N,t,num_times)
m = 1.6*ones(1,d); s = 0.25;
if d > 1
    s = (s^2)*eye(d);
end
if d == 1
    p = @(x) normpdf(log(x), m, s) ./ prod(x,2); 
else
    p = @(x) mvnpdf(log(x), m, s) ./ prod(x,2);
end
lvl = 3:5; %levels 3 through 5
L = 5; % highest level
rho = 0.9; % technical term 1
dt = 0.01; % technical term 2
ellipfn = @(x,l) f(x,l,d); % model operator at highest level, current dim

estimates = zeros(1,num_times);
runtimes = zeros(1,num_times);
for i=1:num_times
   Z = exp(mvnrnd(m,s,N)); %initial log-normal sample 
   tic; 
   [estimates(i),~,~] = MFCE_LogNormal(Z, t, p, ellipfn, lvl, N, rho, dt); 
   runtimes(i) = toc;
end
coeffofvar = mean((estimates - truth).^2) / (truth^2);
avg_time = mean(runtimes);
end