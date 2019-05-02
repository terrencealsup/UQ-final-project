% precompute operators for block with 16

for level = 3:3
    xi = zeros(1,2);
    ACell = {};
    xi(1) = 1;
    [ACell{1}, f] = genOperators2D(@(X, Y)coeffFun2DBlocks2(xi,X, Y), level);
    for i=2:length(xi)
        disp(num2str(i));
        xi = zeros(1, 2);
        xi(i) = 1;
        ACell{i} = genOperators2D(@(X, Y)coeffFun2DBlocks2(xi, X, Y), level);
    end
    
    eval(['save operatorsBlocks2vert_level', num2str(level), '.mat ACell f']);
    
end