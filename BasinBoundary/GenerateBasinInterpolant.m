function [F] = GenerateBasinInterpolant(M)
%GENERATEBASINBOUNDARY 

[Basins] = LoadBasins(M);

if isempty(Basins)
   F = griddedInterpolant; %empty interpolant
   return
end

x = Basins.x(:,1);
y = Basins.x(:,2);

B = ones(size(Basins.B)); %initialize all ones

jj = find(Basins.B~=1); %indeces of basin numbers that are not the initial basin
B(jj) = 3*ones(size(jj)); %all non initial basins are now the number 3

z = B;
A = Basins.A;

L =sqrt(length(x));


xGrid = reshape(x,L,L);
yGrid = reshape(y,L,L);
zGrid = reshape(z,L,L);


F = griddedInterpolant(xGrid,yGrid,zGrid);




end

