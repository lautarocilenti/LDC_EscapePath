function [] = PL_PathEnergy(S,M)
%PL_PATHENERGY 

% plot(M.theta,W)
% thetaSet = M.theta();
% if strcmp(M.rhsString,'Unforced') %bistable case
%     A = [0 1; 3*M.Mrhs.a1 -M.Mrhs.nu]; %a1 already negative so different than write up
% elseif strcmp(M.rhsString,'Duffing') %hardeninng forced case
%     A = [0 1; -3*M.Mrhs.a3*M.qo(1)^2-M.Mrhs.a1, -M.Mrhs.nu]; %a1 already negative so different than write up
% end
% 
% E = [0 0;0 1];
% B = [A E;zeros(size(E)) -A'];
% [ev,e] = eig(B);
% lambda = diag(e);
% iUnstableLambda = find(real(lambda)>0);
% u_ev = ev(:,iUnstableLambda(1)); %unstable eigenvector, just using one of the complex conjugate pair
% xi_q = real(u_ev(1:2)); xi_p = real(u_ev(3:4));
% eta_q = imag(u_ev(1:2)); eta_p = imag(u_ev(3:4));
% T = [xi_q eta_q]\[xi_p eta_p]; %p_eps = T*q_eps
% q_eps = [M.rIC*cos(thetaSet);M.rIC*sin(thetaSet)];
% p_eps = T*q_eps;
% xoSet = [ M.qo.*ones(size(M.qo,1),length(thetaSet))+q_eps; p_eps];
% 
% 
% x = q_eps(1,:);
% y = q_eps(2,:);
% plot3(x,y,S)
% hold on
% [minS,iS] = min(S);
% plot3(x(iS),y(iS),S(iS),'.','markersize',30)
% xlabel('q eps 1')
% ylabel('q eps 2')
% zlabel('S')
% legend('S','minS')
% [xs,is] = sort(x,'ascend');
plot( M.theta,S)
xlabel('$\theta$ (radians)','interpreter','latex')
ylabel('$S(q)$','interpreter','latex')
title(sprintf("Path Energy with $\\nu$: %.2f, R: %.2e",M.Mrhs.nu,M.rIC),'interpreter','latex')

end

