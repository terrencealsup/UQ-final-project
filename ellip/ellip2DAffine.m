function Q = ellip2DAffine(xi, ACell, f, level)

n = size(xi, 1);
Q = zeros(n, 1);
parfor i=1:n
    % assemble system
    A = xi(i, 1)*ACell{1};
    for j=2:length(ACell)
        A = A + xi(i, j)*ACell{j};
    end
    u = A\f;
    
%   p = reshape(u, sqrt(size(u, 1)), sqrt(size(u, 1)));
    Q(i) = getOutput2D(level, reshape(u, sqrt(size(u, 1)), sqrt(size(u, 1))));
end
end