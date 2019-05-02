function p = ellip2DAffine_heat(xi, ACell, f, level)

n = size(xi, 1);
p = zeros(n, 1);
% assemble system
A = xi(:, 1)*ACell{1};
for j=2:length(ACell)
    A = A + xi(:, j)*ACell{j};
end
u = A\f;


p = reshape(u, sqrt(size(u, 1)), sqrt(size(u, 1)));
x = linspace(0,1,size(p,1)); y = linspace(0,1,size(p,2));
[x,y] = meshgrid(x,y);
colormap spring;
contourf(x,y,p,20); colorbar;
xlabel('x'); ylabel('y'); 
end