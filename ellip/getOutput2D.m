function y = getOutput2D(level, u)

assert(size(u, 1) > 1 && size(u, 2) > 1); % make sure we got a matrix; not vector
% integrate over base
h = 1/2^level;
y = h*(u(1, end) + u(end, end)) + 2*h*sum(u(2:end-1, end));

end