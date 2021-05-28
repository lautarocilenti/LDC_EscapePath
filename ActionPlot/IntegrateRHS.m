function [phiSet] = IntegrateRHS(xoSet,M)
%IntegrateRHS - Uses ode45 to integrate Hamiltonian  from initial conditions until
%maximum time or until TerminateAtBoundaryEvent has a value of zero. 

opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1E-2,'Events', M.TerminateEvent);


nt = length(M.tspan(1):M.dt:M.tspan(2));
phiSet = {};
tSet =  zeros(nt,M.nIC); 
parConstant = parallel.pool.Constant(M);
[isCluster] = ProgressBar(size(xoSet,2),"Rise and Fall");
% parfor(ixo = 1:size(xoSet,2),M.nWorkers)
for ixo = 1:size(xoSet,2)
    
    m = parConstant.Value;
    terminateFlag1 = true;
    xo = xoSet(:,ixo);
    m.Mrhs.temp.tOldPeriod = 0;
    m.Mrhs.temp.T = 2*pi/m.Mrhs.w;
    [t,y,~,~,ie] =  m.solver(m.HamiltonianRHS, m.tspan, xo,[opts],m.Mrhs);
    phiSet{ixo} = {t,y,ie};
    if (max(t) >= max(m.tspan))
        error("Integration may have terminated prior to reaching boundary\n")
    end
    if ~isCluster
        fprintf("\b|\n")
    end

    
end
clear parConstant

end



