function [] = PL_MPEP(theta,minPhi,minS,M)
%PL_MPEP 
phi = minPhi{1}{2};
plot(phi(:,1),phi(:,2),'linewidth',2)
hold on
title(sprintf("MPEP with $\\nu$: %.2f, $\\theta$: %.3f, S: %.3f R: %.3e",M.Mrhs.nu,theta,minS,M.rIC),'interpreter','latex')
xlabel("$q_1$",'interpreter','latex')
ylabel("$q_2$",'interpreter','latex')

end

