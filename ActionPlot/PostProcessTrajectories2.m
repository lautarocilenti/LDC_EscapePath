function [phiSetOut] = PostProcessTrajectories2(phiSet,M)
%POSTPROCESSTRAJECTORIES 
if M.plotFall == false
    phiSetOut = phiSet;
    return
end

    phiSetOut = {};

    parConstant = parallel.pool.Constant(M);
    if M.progressbar
        [isCluster] = ProgressBar(length(phiSet),"Post Process Trajectories");
    end
    statusSet = cellfun(@(phi) phi{3},phiSet);
parfor(iPhi = 1:length(phiSet),M.nWorkers)
   m = parConstant.Value;
% for iPhi = 1:length(phiSet)
%     m = M;

       
       t =  phiSet{iPhi}{1}; phi = phiSet{iPhi}{2};
       T = M.Mrhs.T;
       dT = M.dT;
       
       iT = find(mod(t,dT)==0); iT = iT(end-1:end);
       s1Flag = false; s3Flag = false;
       itau = max(iT);
       maxLoopIndex = 1000;
       for loopIndex = 1:maxLoopIndex; %max 50 iterations of loop

           if diff(iT) <= 1
                itau = max(iT);
                tau = t(itau);
                q = phi(itau,1:M.dim)';
                [status,tFall,yFall] = IntegrateToFixedPoint(tau,q,m.Mrhs);
%                 [x1,x2] = InterpFall(m,tFall,yFall);
%                 plot(x1,x2,'o-')
%                 axis([-4,4,-4,4])
                break
           else
               

                tau = t(itau);
                q = phi(itau,1:M.dim)';
                [status,tFall,yFall] = IntegrateToFixedPoint(tau,q,m.Mrhs);
%                 iT
%                 [x1,x2] = InterpFall(m,tFall,yFall);
%                 plot(x1,x2,'o-')
%                 axis([-4,4,-4,4])
 
                if status == m.Mrhs.iA 
%                     if s1Flag
%                         fprintf("Increasing lower bound\n")
                        iT(1) = itau;
%                     end
%                     s1Flag = true;  s3Flag = false;
                elseif status == 2
                    break
                elseif status ~= m.Mrhs.iA  & status ~= 0
%                     if s3Flag
%                         fprintf("Decreasing upper bound\n")
                        iT(end) = itau;
%                     end
%                     s1Flag = false; s3Flag = true;
                else
                    fprintf("Error \n")
                    break
                end
               
           end 
           itau = ceil(mean(iT));

       end
       if loopIndex >= maxLoopIndex
           error("maxLoopIndexReached")
       end
        phiSetOut{iPhi}{1} = t(1:itau,1);
        phiSetOut{iPhi}{2} = phi(1:itau,:);
        phiSetOut{iPhi}{3} = status;
        phiSetOut{iPhi}{4} = tFall;
        phiSetOut{iPhi}{5} = yFall;
        if M.progressbar & ~CheckIfCluster
            fprintf("\b|\n")
        elseif M.progressbar
            fprintf(".")
        end

    end
     fprintf("\n")
    clear parConstant


end


function [x1,x2] = InterpFall(m,tFall,yFall)
    timeinPeriods = 0:2*pi/m.Mrhs.w:max(tFall);
    iTsmall = find(timeinPeriods< min(tFall));
    timeinPeriods = timeinPeriods(iTsmall(end):end);
    x1 = interp1(tFall,yFall(:,1),timeinPeriods);
    x2 = interp1(tFall,yFall(:,2),timeinPeriods);
end
