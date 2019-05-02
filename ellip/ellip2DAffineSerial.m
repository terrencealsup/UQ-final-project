function Q = ellip2DAffineSerial(xi, ACell, f, level)

n = size(xi, 1);
Q = zeros(n, 1);

for i=1:n
    % assemble system
    A = xi(i, 1)*ACell{1};
    for j=2:length(ACell)
        A = A + xi(i, j)*ACell{j};
    end
    u = A\f;
    disp(f);
    
    Q(i) = getOutput2D(level, reshape(u, sqrt(size(u, 1)), sqrt(size(u, 1))));
end