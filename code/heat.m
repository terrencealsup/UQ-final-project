function [Q_samps, xi1, xi2] = heat(dim,n_points,level)
% function handle to take visualize the dependence of f on xi in dimensions
% 1 and 2, for UNIFORM xi, as a plot / contour 
%
% Args: dim -- number of dimensions in parameter space, either 1 or 2
%       n_points -- number of points to put down on the interval [1,10] in
%                   each dimension
%       level -- level for the PDE solution

addpath(genpath('ellip')); %add ellip subdirectory to path

if dim == 1
    xi1 = linspace(0.5,5.5,n_points)';
    xi2 = [];
    xi = xi1;
else
    m = n_points^2;
    
    xi1 = linspace(0.5,5.5,n_points);
    xi2 = linspace(0.5,5.5,n_points);
    [xi1,xi2] = meshgrid(xi1,xi2); %mesh grid
    xi = [reshape(xi1,m,1), reshape(xi2,m,1)]; %reshape for ellip2DAffine
end

ops = load(['operatorsBlocks',num2str(dim),'_level', num2str(level), '.mat']);  %load ops
Q_samps = ellip2DAffine(xi, ops.ACell, ops.f, level);

if dim == 2
   Q_samps = reshape(Q_samps,n_points,n_points); 
end

if dim == 1
    plot(xi1,Q_samps,'-s','LineWidth',2);
    xlabel('$\xi$','Interpreter','latex'); ylabel('QoI'); title('QoI vs $\xi$','Interpreter','latex');
    xlim([0.5,5]); xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5]);
    set(gca,'FontSize',16);
else
%     contourf(xi1,xi2,Q_samps,0.1:0.1:1); colorbar;
    contourf(xi1,xi2,Q_samps); colorbar;
    xlabel('$\xi_{1}$','Interpreter','latex'); ylabel('$\xi_{2}$','Interpreter','latex');
    title('QoI vs $\xi_{1}, \xi_{2}$','Interpreter', 'latex');
    xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5]); yticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5]);
    set(gca,'FontSize',16);
end

end