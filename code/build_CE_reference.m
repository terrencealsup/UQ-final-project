close all; clearvars;
% compute small probability rare events using CE
addpath('../ellip/');

d = 16; %dimension
N = 1e5; %number of samples
t = 0.3; %threshold (see above)
tname = num2str(t); tname(1:2) = []; %can't have the 0. in file name

disp(['Building d=',num2str(d),' high-fidelity CE with N=', num2str(N),...
    ' for threshold t=',num2str(t)]);

%parameters on the target lognormal
m = 1.6*ones(1,d); s = 0.25;
if d > 1
    s = (s^2)*eye(d);
end

%density of the target log normal
if d == 1
    p = @(x) normpdf(log(x), m, s) ./ prod(x,2); 
else
    p = @(x) mvnpdf(log(x), m, s) ./ prod(x,2);
end

Z = exp(mvnrnd(m,s,N)); %initial log-normal sample

lvl = 3:5; %levels 3 through 5
L = 5; % highest level
rho = 0.9; % technical term 1
dt = 0.01; % technical term 2

ellipfn = @(x) f(x,L,d); % model operator at highest level, current dim

tic; [prob, mu, S] = CE_LogNormal(Z, t, p, ellipfn, N, rho, dt,{}); timed = toc;

CE_truth = prob;
disp(['high-fidelity CE with N=', num2str(N), ' samples took t=',num2str(timed),' seconds']);

disp(['CE estimate is ', num2str(prob)]);
str = input('Do you want to save current reference? y / n: ', 's');
% str = 'y'; %autosave
if str == 'y'
    save(['CEdim',num2str(d),'threshold0point',num2str(tname),'.mat'],'CE_truth');
    disp('Saved reference');
else
    disp('Did not save reference');
end





function vals = f(x,l,d)
    ops = load(['operatorsBlocks',num2str(d),'_level',num2str(l),'.mat']);
    vals = ellip2DAffine(x,ops.ACell,ops.f,l);
end