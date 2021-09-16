function [phiSet] = IntegrateRHS(xoSet,M)
%IntegrateRHS - Uses ode45 to integrate Hamiltonian  from initial conditions until
%maximum time or until TerminateAtBoundaryEvent has a value of zero. 

opts = odeset('RelTol',1e-6,'AbsTol',1e-6);


nt = length(M.tspan(1):M.dt:M.tspan(2));
phiSet = {};
tSet =  zeros(nt,M.nIC); 
parConstant = parallel.pool.Constant(M);
if M.progressbar
    [isCluster] = ProgressBar(size(xoSet,2),"Rise and Fall");
end
parfor(ixo = 1:size(xoSet,2),M.nWorkers)
    m = parConstant.Value;
% for ixo = 1:size(xoSet,2)
% m = M;
    status = []; t= [];

    
    
    xo = xoSet(:,ixo);
%     A = m.Mrhs.FixedPoints.FP(m.Mrhs.iA,:);
    T = 2*pi/m.Mrhs.w;
    to = m.pp/m.Mrhs.w;
    dT = m.dT;
    tVector = 0:dT:m.tspan(end);
    tVector(1) = to;
    tAll = [tVector(1)]; yAll = [xo'];
    iT = 1;
    while  iT <=length(tVector)
        tspan = [tVector(iT);tVector(iT+1)];

        [t,y] =  m.solver(m.HamiltonianRHS, tspan, xo,[opts],m.Mrhs);
        tAll = [tAll; t(2:end)];
        yAll = [yAll;y(2:end,:)];
        
        A = m.Mrhs.FixedPoints.GetFixedPoint(mod(t(end),T),m.Mrhs.iA);
        if norm(y(end,1:M.dim)) > 1E3 %solution blew up
            if dT < T/8
                fprintf("Terminating because dT is too small\n")
                break 
            end
            fprintf("Solution blew up, starting over,dT = %f\n",dT/2)
            dT = dT/2;
            tVector = 0:dT:m.tspan(end);
            tVector(1) = to;
            xo = xoSet(:,ixo);
            tAll = [tVector(1)]; yAll = [xo'];
            iT = 1;
            continue
        elseif norm(y(end,1:M.dim)-A)>=m.Mrhs.rA1 
            [status,~,~] = IntegrateToFixedPoint(t(end),y(end,1:M.dim),m.Mrhs);
            
            if status~=m.Mrhs.iA & status~=0
               break 
            elseif status == 0 
                fprintf("Error, status = 0 \n")
            end
        end
        xo = y(end,:)';
        iT = iT+1;
    end
%     status
    phiSet{ixo} = {tAll,yAll,status};
    if (max(t) >= max(m.tspan))
        error("Integration may have terminated prior to reaching boundary\n")
    end
    if M.progressbar & ~CheckIfCluster
        fprintf("\b|\n")
    elseif M.progressbar
        fprintf(".")
    end

    
end
fprintf("\n")
clear parConstant

end



