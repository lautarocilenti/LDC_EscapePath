function [S] = IntegrateLagrangian(phiSet,M)
%INTEGRATELAGRANGIAN 

S = zeros(1,length(phiSet));

for iPhi = 1:length(phiSet)
   phi = phiSet{iPhi}{2};
   t = phiSet{iPhi}{1};
   S(1,iPhi) = trapz(t,M.Lagrangian(phi,M));
end


end

