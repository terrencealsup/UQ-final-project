function Z = coeffFun2DBlocks2(xi, X, Y)
% In
%   xi      ...     parameter
%   X, Y    ...     coordinates (with meshgrid)
% Out
%   coeff   ...     coefficient

k1 = 2;
k2 = 1;
h1 = 1/k1;
h2 = 1/k2;
Z = zeros(size(X));
for i=1:k1
    for j=1:k2
        if(j == k2 && i < k1)
            Z((i-1)*h1 <= X & X < i*h1 & (j-1)*h2 <= Y & Y <= j*h2) = xi(k1*(i - 1) + j);
        elseif(i ==k1 && j < k2)
            Z((i-1)*h1 <= X & X <= i*h1 & (j-1)*h2 <= Y & Y < j*h2) = xi(k1*(i - 1) + j);
        elseif(i ==k1 && j == k2)
            Z((i-1)*h1 <= X & X <= i*h1 & (j-1)*h2 <= Y & Y <= j*h2) = xi(k1*(i - 1) + j);
        else
            Z((i-1)*h1 <= X & X < i*h1 & (j-1)*h2 <= Y & Y < j*h2) = xi(k1*(i - 1) + j);
        end
        
    end
end

end