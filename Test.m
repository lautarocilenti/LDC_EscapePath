function [] = Test()

M.a1 = 1;
M.a3 = .3; 
M.nu = .1;
M.F = .4; 
M.w = 1.4; 
M.kc = .01; 
M.wc = .1;

%2*pi;%rhs parameters
qo = [1.3586;2.4173];
M.filterDim = 1;
Mrhs = M;
Xj = [qo',zeros(1,M.filterDim),zeros(size(qo')),zeros(1,M.filterDim)];
J =DuffingFilteredJacobian(Xj,Mrhs);

T = 2*pi/M.w;
xo = [0;0;0;0;0;1];
[t,y] = ode45(@(t,x) J*x,[0 T],xo)

plot(t,y)

end

function [J] = DuffingFilteredJacobian(qo,Mrhs)
        A = [0 1 0;...
            -3*Mrhs.a3*qo(1)^2-Mrhs.a1, -Mrhs.nu 1;
            0 0 -Mrhs.wc]; %linearized original system at attractor coordinate
        E = [0 0 0;0 0 0; 0 0 Mrhs.wc^2]; %contribution of p to original system
        J = [A E;zeros(size(E)) -A.']; %Jacobian of hamiltonian system 
end

