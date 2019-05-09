close all; clearvars;
% Testing MFCE on our elliptic PDE model problem
addpath('../ellip/');
%test against MC truths
%later test against the CE truths for smaller rare events
% d = 1: -- t = 0.40 --> Pt = 3.13e-03
%        -- t = 0.46 --> Pt = 5.12e-04
%        -- t = 0.50 --> Pt = 1.45e-04
%  
% d = 4: -- t = 0.30 --> Pt = 3.30e-03
%        -- t = 0.32 --> Pt = 8.28e-04
%        -- t = 0.34 --> Pt = 1.83e-04
% 
% d = 16: -- t = 0.255 --> Pt = 1.05e-03
%         -- t = 0.260 --> Pt = 4.48e-04
%         -- t = 0.265 --> Pt = 1.75e-04

d = 16; %dimension
N = 1e4; %number of samples
t = 0.265; %threshold (see above)
%change truth based on above
truth = 1.75e-04;

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

ellipfn = @(x,l) f(x,l,d); % model operator at highest level, current dim

tic; [prob, mu, S] = MFCE_LogNormal(Z, t, p, ellipfn, lvl, N, rho, dt); timed = toc;

rel = abs(prob-truth)/truth;

disp(['CE estimate is ', num2str(prob)]);
disp(['truth is ',num2str(truth)]);
disp(['relative error is ', num2str(rel)]);
disp(['High fidelity CE using N=',num2str(N), ' took ',num2str(timed), ' seconds']);



function vals = f(x,l,d)
    ops = load(['operatorsBlocks',num2str(d),'_level',num2str(l),'.mat']);
    vals = ellip2DAffine(x,ops.ACell,ops.f,l);
end