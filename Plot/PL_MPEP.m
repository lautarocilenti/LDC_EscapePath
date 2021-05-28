function [iS] = PL_MPEP(data)
%PL_MPEP 
theta = data.theta(data.minPhiIndex);
minPhi = data.phiSet{data.minPhiIndex}; 
minS = data.minS;
M = data.M;



[minS,iS] = min(minS);
t = minPhi{1};
phi = minPhi{2};
tFall = minPhi{4};
phiFall = minPhi{5};
if strcmp(M.rhsString,'Unforced');
    plot(phi(:,1),phi(:,2),'linewidth',2)
else strcmp(M.rhsString,'Duffing');
    color = [rand() rand() rand()];
    it = find(mod(M.Mrhs.w*t,2*pi)<M.Mrhs.psiEps);
    plot(phi(it,1),phi(it,2),'o-','Color',color);
    if ~strcmp(M.terminateType,"Saddle")
        itF = find(mod(M.Mrhs.w*tFall,2*pi)<M.Mrhs.psiEps);
        plot([phi(it(end),1);phiFall(itF,1)],[phi(it(end),2);phiFall(itF,2)],'x--','Color',color);
    end
end
hold on
title(sprintf("MPEP with $\\nu$: %.2f, $\\theta$: %.3f, S: %.3f R: %.3e",M.Mrhs.nu,theta,minS,M.rIC),'interpreter','latex')
xlabel("$q_1$",'interpreter','latex')
ylabel("$q_2$",'interpreter','latex')

end

