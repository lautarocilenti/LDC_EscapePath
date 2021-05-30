function [] = PL_PathEnergy(data)
%PL_PATHENERGY 


plot( data.theta,data.S)
xlabel('$\theta$ (radians)','interpreter','latex')
ylabel('$S(q)$','interpreter','latex')
title(sprintf("Path Energy with $\\nu$: %.2f, R: %.2e",data.M.Mrhs.nu,data.M.rIC),'interpreter','latex')

end

