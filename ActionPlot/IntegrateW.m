function [W] = IntegrateW(tSet,phiSet,iIntersect,M)
%INTEGRATELAGRANGIAN 

W = zeros(M.nIC,1);

for iPhi = 1:M.nIC
   phi = phiSet(:,:,iPhi);
   iTf = iIntersect(iPhi);
    if iTf == -1
        fprintf("Path with initial ccndition [%f,%f], did not intersect boundary\n",phi(1,1),phi(1,2));
        W(iPhi) = -1;
    else
        W(iPhi) = trapz(tSet(1:iTf,iPhi),Wd(phi(1:iTf,3:4)));
    end
end
end

function [rhs] = Wd(p)
    rhs = 1/2*(p(:,1).^2+p(:,2).^2);
end