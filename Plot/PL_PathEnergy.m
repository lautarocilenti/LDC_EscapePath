function [] = PL_PathEnergy(data)
%PL_PATHENERGY 
data

% phiSet = data.phiSet;
% for i = 1:length(phiSet)
%    p(i,:) = phiSet{i}{2}(1,1:2)-data.attractors(1,:); 
% end
% plot( p(:,1), p(:,2),'s-r')
[theta, isort] = sort(mod(data.theta,2*pi),'ascend');
plot( theta,data.S(isort))
xlabel('$\theta$ (radians)','interpreter','latex')
ylabel('$S$','interpreter','latex')
%title(sprintf("Path Energy with $\\nu$: %.2f, R: %.2e",data.M.Mrhs.nu,data.M.rIC),'interpreter','latex')

end

