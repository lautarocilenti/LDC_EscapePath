function [xoSet] = GenerateInitialConditions(thetaSet,M)
%GENERATEINITIALCONDITIONS Creates set of initial conditions
%
% Input: thetaSet: set of angles surrounding fixed point of poincare section
%        M: large  structure of parameters
%
% Output: xoSet: each column contains 2D values and is an initial condition
%                of the 2D hamiltonian system


%objects for rhs generation 
    x = sym('x',[M.dim*2,1]);
    syms tau
    
%Grab the location of the attractor in the poincare section
[qo] = GetInitialFixedPoint(M.pp,M.Mrhs);


if strcmp(M.rhsString,'Duffing') %Duffing forced case
    A = [0 1; -3*M.Mrhs.a3*x(1)^2-M.Mrhs.a1, -M.Mrhs.nu]; %linearized original system
    E = [0 0;0 1]; %contribution of p to original system
    J = [A E;zeros(size(E)) -A.']; %Jacobian of hamiltonian system
else
    error("Unknown Jacobian\n")
end

jRHS = matlabFunction(J*x,'vars',{tau,x}); %convert to function handle for integration

tspan = [0,M.Mrhs.T]; %one period time span
I = eye(M.dim*2); %identity matrix
opts = odeset('RelTol',1e-12,'AbsTol',1e-12); %ode45 options
Psi_T = zeros(size(I)); %initialize Psi

for i = 1:M.dim*2
    xo = I(:,i); %initial condition is identity matrix columns
    [t,y] =  M.solver(jRHS, tspan, xo,[opts]);
    Psi_T(:,i) = y(end,:).'; %columns of Psi at time T
end

E = Psi_T;

%get eigenvalues
[ev,e] = eig(E);
lambda = diag(e); %eigenvalues 
[lambda,isortlambda] = sort(vecnorm(lambda,2,2),'ascend'); %sort eigenvalues by distance from origin in real imag plane
iUnstableLambda = find(lambda>1); 
iStableLambda = find(lambda<=1);

ev = ev(:,isortlambda); %re sort eigenvectors 

u_ev = ev(:,iUnstableLambda(1)); %eigenvectors of unstable eigenvalues, just grab 1 of conjugate pair

Z_xv = [real(u_ev(1:M.dim,:)),imag(u_ev(1:M.dim,:))];
Z_pxpv = [real(u_ev(M.dim+1:2*M.dim,:)),imag(u_ev(M.dim+1:2*M.dim,:))];
T = Z_xv\Z_pxpv; %transformation matrix


q_eps = [M.rIC*cos(thetaSet);M.rIC*sin(thetaSet)]; %initial offset from attractor

p_eps = T*q_eps; %initial momenta

xoSet = [qo.*ones(size(qo,1),length(thetaSet))+q_eps;p_eps];





end

function [qo] = GetInitialFixedPoint(pp,Mrhs)
    t = pp/Mrhs.w;
    [fp,~] = Mrhs.FixedPoints.GetFixedPoint(mod(t,Mrhs.T),Mrhs.iA);
    qo = fp';
end