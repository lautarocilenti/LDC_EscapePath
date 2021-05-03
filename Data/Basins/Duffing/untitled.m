clear all
filename = "Basins_Phi00.mat";
load(filename);

Basins.A
iList = [];
[imin,iminIndex] = min(iList);
iList(iminIndex) = [];

for j = 1:length(iList)
    i = iList(j);
    ii = find(Basins.B == i);
    Basins.B(ii) = imin;
end

k = size(Basins.A,1); 
kList = 1:k;
k2 = ismember(kList,iList);
A = Basins.A(~k2,:);

Basins.A = A;
Basins.iH = 2;
Basins.A
save(filename,'Basins');