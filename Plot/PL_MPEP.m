function [iS] = PL_MPEP(data)
%PL_MPEP 
theta = data.minTheta
minPhi = data.minPhi; 
minS = data.minS;
M = data.M;
psi = data.psi;


[minS,iS] = min(minS);
t = minPhi{iS}{1};
phi = minPhi{iS}{2};
if strcmp(M.rhsString,'Unforced')
    plot(phi(:,1),phi(:,2),'linewidth',2)
else strcmp(M.rhsString,'Duffing')
    it = find(abs(psi{iS}-mod(M.Mrhs.w*t,2*pi))<M.Mrhs.psiEps);
    plot(phi(it,1),phi(it,2),'linewidth',2)
end
hold on
title(sprintf("MPEP with $\\nu$: %.2f, $\\theta$: %.3f, S: %.3f R: %.3e",M.Mrhs.nu,theta,minS,M.rIC),'interpreter','latex')
xlabel("$q_1$",'interpreter','latex')
ylabel("$q_2$",'interpreter','latex')

end

