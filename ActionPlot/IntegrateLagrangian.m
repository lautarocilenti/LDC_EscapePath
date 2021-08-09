function [S,pNorm] = IntegrateLagrangian(phiSet,M)
%INTEGRATELAGRANGIAN 

S = zeros(length(phiSet),1);

for iPhi = 1:length(phiSet)
   phi = phiSet{iPhi}{2};
   pNorm(iPhi) = norm(phi(end,3:4));
   t = phiSet{iPhi}{1};
   S(iPhi) = trapz(t,M.Lagrangian(phi,M));
end


end

