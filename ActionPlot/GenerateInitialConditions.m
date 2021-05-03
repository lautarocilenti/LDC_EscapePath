function [xoSet] = GenerateInitialConditions(thetaSet,M)
%GENERATEINITIALCONDITIONS Creates set of initial momenta for integration
%of the Hamiltonian space.


if strcmp(M.rhsString,'Unforced') %bistable case
    A = [0 1; 3*M.Mrhs.a1 -M.Mrhs.nu]; %a1 already negative so different than write up
elseif strcmp(M.rhsString,'Duffing') %hardeninng forced case
    A = [0 1; -3*M.Mrhs.a3*M.qo(1)^2-M.Mrhs.a1, -M.Mrhs.nu]; %a1 already negative so different than write up
end

E = [0 0;0 1];
B = [A E;zeros(size(E)) -A'];
[ev,e] = eig(B);
lambda = diag(e);
iUnstableLambda = find(real(lambda)>0);
u_ev = ev(:,iUnstableLambda(1)); %unstable eigenvector, just using one of the complex conjugate pair
xi_q = real(u_ev(1:2)); xi_p = real(u_ev(3:4));
eta_q = imag(u_ev(1:2)); eta_p = imag(u_ev(3:4));
T = [xi_q eta_q]\[xi_p eta_p]; %p_eps = T*q_eps
q_eps = [M.rIC*cos(thetaSet);M.rIC*sin(thetaSet)];
p_eps = T*q_eps;
xoSet = [ M.qo.*ones(size(M.qo,1),length(thetaSet))+q_eps; p_eps];




end

