function [] = PL_PathEnergy(data)
%PL_PATHENERGY 
data

phiSet = data.phiSet;
for i = 1:length(phiSet)
   p(i,:) = phiSet{i}{2}(1,1:2)-data.attractors(1,:); 
end
plot( p(:,1), p(:,2),'s-r')
% plot( data.theta,data.S)
xlabel('$\theta$ (radians)','interpreter','latex')
ylabel('$S(q)$','interpreter','latex')
%title(sprintf("Path Energy with $\\nu$: %.2f, R: %.2e",data.M.Mrhs.nu,data.M.rIC),'interpreter','latex')

end

