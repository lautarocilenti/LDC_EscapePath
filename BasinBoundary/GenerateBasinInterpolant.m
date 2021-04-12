function [F,A] = GenerateBasinInterpolant(M)
%GENERATEBASINBOUNDARY 

[Basins] = LoadBasins(M);

x = Basins.x(:,1);
y = Basins.x(:,2);

ii = find(Basins.B~=M.BN); %indeces of basin numbers that are not the initial basin
Basins.B(ii) = 3*ones(size(ii)); %all non initial basins are now the number 3

z = Basins.B;
A = Basins.A;

L =sqrt(length(x));


xGrid = reshape(x,L,L);
yGrid = reshape(y,L,L);
zGrid = reshape(z,L,L);

F = griddedInterpolant(xGrid,yGrid,zGrid);





end

