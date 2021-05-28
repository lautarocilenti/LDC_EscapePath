function [phiSetOut] = PostProcessTrajectories(phiSet,M)
%POSTPROCESSTRAJECTORIES 
phiSetOut = {};
if strcmp(M.terminateType,"Saddle")
    parConstant = parallel.pool.Constant(M);
    [isCluster] = ProgressBar(length(phiSet),"Post Process Trajectories");
    status = cellfun(@(phi) phi{3},phiSet)
    phiSet = phiSet(find(status)==1);
%    
%     parfor(iPhi = 1:length(phiSet),Inf)
    for iPhi = 1:length(phiSet)
       m = parConstant.Value;
       t =  phiSet{iPhi}{1}; phi = phiSet{iPhi}{2};
               for i = length(t):-1:1
                  itau = i;
                  tau = t(itau);
                  q = phi(itau,1:2)';
                  [status,tFall,yFall] = IntegrateToFixedPoint(tau,q,m.Mrhs);
                  if status == 1
                      if itau < length(t)
                      itau = itau+1;
                      end
                      phiSetOut{iPhi}{1} = t(1:itau,1);
                      phiSetOut{iPhi}{2} = phi(1:itau,:);
                      phiSetOut{iPhi}{3} = phiSet{iPhi}{3};
                      tau = t(itau);
                      q = phi(itau,1:2)';
                      [status,tFall,yFall] = IntegrateToFixedPoint(tau,q,m.Mrhs);
                      phiSetOut{iPhi}{3} = status;
                      phiSetOut{iPhi}{4} = tFall;
                      phiSetOut{iPhi}{5} = yFall;
                      break
                  end

               end

           
      if ~isCluster
            fprintf("\b|\n")
       end

    end
    clear parConstant
else
    parConstant = parallel.pool.Constant(M);
    [isCluster] = ProgressBar(length(phiSet),"Post Process Trajectories");
    parfor(iPhi = 1:length(phiSet),Inf)
       m = parConstant.Value;
       t =  phiSet{iPhi}{1}; phi = phiSet{iPhi}{2};
       for i = length(t):-1:1
          itau = i;
          tau = t(itau);
          q = phi(itau,1:2)';
          [status,tFall,yFall] = IntegrateToFixedPoint(tau,q,m.Mrhs);
          if status == 1
              itau = itau+1;
              phiSetOut{iPhi}{1} = t(1:itau,1);
              phiSetOut{iPhi}{2} = phi(1:itau,:);
              phiSetOut{iPhi}{3} = phiSet{iPhi}{3};
              tau = t(itau);
              q = phi(itau,1:2)';
              [status,tFall,yFall] = IntegrateToFixedPoint(tau,q,m.Mrhs);
              phiSetOut{iPhi}{3} = status;
              phiSetOut{iPhi}{4} = tFall;
              phiSetOut{iPhi}{5} = yFall;
              break
          end

       end
       if ~isCluster
            fprintf("\b|\n")
       end

    end
    clear parConstant
end

end

