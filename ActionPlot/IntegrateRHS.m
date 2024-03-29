function [phiSet] = IntegrateRHS(xoSet,M)
%IntegrateRHS - Uses ode45 to integrate Hamiltonian  from initial conditions until
%maximum time or until TerminateAtBoundaryEvent has a value of zero. 
warning off
if M.includePhase
   phaseSet = xoSet(end,:);
   xoSet(end,:) = [];
else
    phaseSet = M.pp*ones(size(xoSet(end,:)));
end


opts = odeset('RelTol',1e-6,'AbsTol',1e-6);


nt = length(M.tspan(1):M.dt:M.tspan(2));
phiSet = {};
tSet =  zeros(nt,M.nIC); 
parConstant = parallel.pool.Constant(M);
if M.progressbar
    [isCluster] = ProgressBar(size(xoSet,2),"Rise and Fall");
end

% 
parfor(ixo = 1:size(xoSet,2))
    m = parConstant.Value;
%     for ixo = 1:size(xoSet,2)
%     m = M;
        status = []; t= [];



        xo = xoSet(:,ixo);
    %     A = m.Mrhs.FixedPoints.FP(m.Mrhs.iA,:);
        T = 2*pi/m.Mrhs.w;
        if m.includePhase
           to = phaseSet(ixo)/m.Mrhs.w;
        else
           to = m.pp/m.Mrhs.w;
        end
        
        dT = m.dT;
        tVector = 0:dT:m.tspan(end);
        iSmaller = tVector <= to;
        tVector(iSmaller) = [];
        tVector = [to tVector];
        tAll = [tVector(1)]; yAll = [xo'];
        iT = 1;
        while  iT <length(tVector)
            tspan = [tVector(iT);tVector(iT+1)];

            [t,y] =  m.solver(m.HamiltonianRHS, tspan, xo,[opts],m.Mrhs);
            tAll = [tAll; t(2:end)];
            yAll = [yAll;y(2:end,:)];

            A = m.Mrhs.FixedPoints.GetFixedPoint(mod(t(end),T),m.Mrhs.iA);
    
            if m.dim == 2
                blewUpThreshold = 1E2;
            elseif M.dim <=6
                blewUpThreshold = 1E2;
            elseif M.dim > 6
                blewUpThreshold = 1E2;
            end


    %         if norm(y(end,1:M.dim)) > blewUpThreshold %solution blew up
            if norm(y(end,:)) > blewUpThreshold 
               if m.Mrhs.a3 < 0 
                   status = size(m.Mrhs.FixedPoints.FP,1)+1;%infinity has a basin
                   break
               end
               if dT < T/4
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
                    errID = 'MyComponent:noSuchVariable';
                    msg = ['StatusIs0,Position:',sprintf(repmat('%.4f,',1,M.dim),y(end,1:M.dim))];
                    baseException = MException(errID,msg);
%                     throw(baseException)
                end
            end
            xo = y(end,:)';
            iT = iT+1;
        end
    %     status
        if isempty(status)
            dummy = 1;
        end
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



