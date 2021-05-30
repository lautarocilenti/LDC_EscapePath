function [phiSetOut] = PostProcessTrajectories2(phiSet,M)
%POSTPROCESSTRAJECTORIES 
if M.plotFall == false
    phiSetOut = phiSet;
    return
end

    phiSetOut = {};

    parConstant = parallel.pool.Constant(M);
    [isCluster] = ProgressBar(length(phiSet),"Post Process Trajectories");
    statusSet = cellfun(@(phi) phi{3},phiSet);
    parfor(iPhi = 1:length(phiSet),M.nWorkers)
% for iPhi = 1:length(phiSet)
       m = parConstant.Value;
       t =  phiSet{iPhi}{1}; phi = phiSet{iPhi}{2};
       T = M.Mrhs.T;
       iT = find(mod(t,T)==0); iT = iT(end-1:end);
      
       for loopIndex = 1:50; %max 50 iterations of loop
           
           if diff(iT) <= 1
                break
           else
               itau = ceil(mean(iT));
                tau = t(itau);
                q = phi(itau,1:2)';
                [status,tFall,yFall] = IntegrateToFixedPoint(tau,q,m.Mrhs);
                
                if status == 1
                    iT(1) = itau;
                elseif status == 2
                    break
                else
                    iT(end) = itau;
                end
               
           end 
       end
        phiSetOut{iPhi}{1} = t(1:itau,1);
        phiSetOut{iPhi}{2} = phi(1:itau,:);
        phiSetOut{iPhi}{3} = status;
        phiSetOut{iPhi}{4} = tFall;
        phiSetOut{iPhi}{5} = yFall;
       if ~isCluster
            fprintf("\b|\n")
       end

    end
    clear parConstant


end

