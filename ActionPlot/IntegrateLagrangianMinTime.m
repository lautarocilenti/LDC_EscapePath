function [S] = IntegrateLagrangianMinTime(phiSet,tmax,M)
%INTEGRATELAGRANGIAN 

S = zeros(length(phiSet),1);

for iPhi = 1:length(phiSet)
   phi = phiSet{iPhi}{2};
   t = phiSet{iPhi}{1};
   
   if t(end) <= tmax
        S(iPhi) = trapz(t,M.Lagrangian(phi,M));
   else
       ii = find(t<tmax);
       for d = 1:size(phi,2)
            phiEnd(1,d) = interp1(t,phi(:,d),tmax);
       end
       
       tNew = [t(ii);tmax];
       phiNew = [phi(ii,:);phiEnd];
       S(iPhi) = trapz(tNew,M.Lagrangian(phiNew,M));
   end
end


end

