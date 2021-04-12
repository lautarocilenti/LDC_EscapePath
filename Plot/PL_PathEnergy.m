function [] = PL_PathEnergy(S,M)
%PL_PATHENERGY 

% plot(M.theta,W)
semilogy(M.theta,S)
xlabel('$\theta$ (radians)','interpreter','latex')
ylabel('$S(q)$','interpreter','latex')
title(sprintf("Path Energy with $\\nu$: %.2f, R: %.2e",M.Mrhs.nu,M.rIC),'interpreter','latex')

end

