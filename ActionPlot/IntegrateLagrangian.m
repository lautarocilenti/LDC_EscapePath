function [S,phiSet] = IntegrateLagrangian(phiSet,M)
%INTEGRATELAGRANGIAN 

S = zeros(1,length(phiSet));

for iPhi = 1:length(phiSet)
   phi = phiSet{iPhi}{2};
   fA = phiSet{iPhi}{3};
   t = phiSet{iPhi}{1};
   S(1,iPhi) = trapz(t,M.Lagrangian(phi,M));
   if M.Mrhs.fA ~= 0 & fA ~= M.Mrhs.fA
       S(1,iPhi) = S(1,iPhi)+10;
   end
   phiSet{iPhi}{6} = S(1,iPhi);
end


end

