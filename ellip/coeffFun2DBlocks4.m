function Z = coeffFun2DBlocks4(xi, X, Y)
% In
%   xi      ...     parameter
%   X, Y    ...     coordinates (with meshgrid)
% Out
%   coeff   ...     coefficient

k = 2;
h = 1/k;
Z = zeros(size(X));
for i=1:k
    for j=1:k
        if(j == k && i < k)
            Z((i-1)*h <= X & X < i*h & (j-1)*h <= Y & Y <= j*h) = xi(k*(i - 1) + j);
        elseif(i ==k && j < k)
            Z((i-1)*h <= X & X <= i*h & (j-1)*h <= Y & Y < j*h) = xi(k*(i - 1) + j);
        elseif(i ==k && j == k)
            Z((i-1)*h <= X & X <= i*h & (j-1)*h <= Y & Y <= j*h) = xi(k*(i - 1) + j);
        else
            Z((i-1)*h <= X & X < i*h & (j-1)*h <= Y & Y < j*h) = xi(k*(i - 1) + j);
        end
        
    end
end

end