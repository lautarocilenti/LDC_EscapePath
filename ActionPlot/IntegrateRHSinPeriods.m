function [phiSet] = IntegrateRHS(xoSet,M)
%IntegrateRHS - Uses ode45 to integrate Hamiltonian  from initial conditions until
%maximum time or until TerminateAtBoundaryEvent has a value of zero. 

opts = odeset('RelTol',1e-6,'AbsTol',1e-6);


nt = length(M.tspan(1):M.dt:M.tspan(2));
phiSet = {};
tSet =  zeros(nt,M.nIC); 
parConstant = parallel.pool.Constant(M);
[isCluster] = ProgressBar(size(xoSet,2),"Rise and Fall");
parfor(ixo = 1:size(xoSet,2),M.nWorkers)
% for ixo = 1:size(xoSet,2)
    
    m = parConstant.Value;
    terminateFlag1 = true;
    xo = xoSet(:,ixo);
    T = 2*pi/m.Mrhs.w;
    tVector = 0:T:m.tspan(end);
    tAll = [tVector(1)]; yAll = [xo'];
    for iT = 1:length(tVector)
        tspan = [tVector(iT);tVector(iT+1)];

        [t,y] =  m.solver(m.HamiltonianRHS, tspan, xo,[opts],m.Mrhs);
        tAll = [tAll; t(2:end)];
        yAll = [yAll;y(2:end,:)];
        if norm(y(end,1:M.dim)-m.Mrhs.qo')>=m.Mrhs.rA1 
            [status,t1,y1] = IntegrateToFixedPoint(t(end),y(end,1:M.dim),m.Mrhs);

            if status>1
               break 
            end
        end
        xo = y(end,:)';
    end
    phiSet{ixo} = {tAll,yAll,status};
    if (max(t) >= max(m.tspan))
        error("Integration may have terminated prior to reaching boundary\n")
    end
    if ~isCluster
        fprintf("\b|\n")
    end

    
end
clear parConstant

end



