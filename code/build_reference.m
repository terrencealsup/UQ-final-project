clear variables; close all;
%computes references using MC
addpath('../ellip/');


m = 1e6; %number of samples for MC
% num_blocks = 8; %number of blocks
l = 5; %reference level

for num_blocks=[16]
disp(['Generating m=', num2str(m), ' level ',num2str(l), ' samples with ',...
    'd=',num2str(num_blocks),' many blocks']);


% pd = makedist('Lognormal', 'mu', 1.6,'sigma',0.25);
% xi = random(pd,m,1);
mu = 1.6;
sigma = 0.25;
xi = exp(mvnrnd(mu*ones(1,num_blocks), (sigma^2)*eye(num_blocks),m));
ops = load(['operatorsBlocks',num2str(num_blocks),'_level', num2str(l), '.mat']);  %load ops
tic; Q_samps = ellip2DAffine(xi, ops.ACell, ops.f, l); t = toc;%samples Q m times
disp(['m=', num2str(m), ' level ',num2str(l), ' samples with ',...
    'd=',num2str(num_blocks),' many blocks took t=',num2str(t),' seconds']);

str = input('Do you want to save current samples? y / n: ', 's');
% str = 'y'; %autosave
if str == 'y'
    save(['LogNormSampsBlocks',num2str(num_blocks),'_level', num2str(l), '.mat'],'Q_samps');
    disp('Saved samples');
else
    disp('Did not save samples');
end
clear Q_samps
end






function Q_mean = Vanilla_MC(m,num_blocks,l)
% Computes an MC estimate of Q(u) = int_(base) u(x) dx
% using the level l, with m samples
% 
% Args: -- m : the number samples
%       -- l : the level to sample at
% 

%load operators at level l -> corresponds to PDE solve and quadrature
%weights for the QoI at base.
ops = load(['operatorsBlocks',num2str(num_blocks),'_level', num2str(l), '.mat']); 

xi = rand(m, 16)*(10 - 1) + 1; %sample the RF for conductivity m times
Q_samps = ellip2DAffine(xi, ops.ACell, ops.f, l); %samples Q m times
Q_mean = mean(Q_samps);
end