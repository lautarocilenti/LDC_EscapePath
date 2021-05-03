function [phiSet2,bIOut,psiOut] = PostProcessTrajectories(phiSet,M)
%POSTPROCESSTRAJECTORIES 
if strcmp(M.rhsString,"Duffing")
    biSet = M.biSet{1};
    phaseSet = M.biSet{2};

    for iPhi = 1:length(phiSet)
       phi = phiSet{iPhi}{2};
       t = phiSet{iPhi}{1};
       exitFlag = 0;
       L = length(t);
       for it = L:-1:1
           tau = t(it);
           y = phi(it,:);
           psi = M.Mrhs.w*tau;
           for i = length(phaseSet)
               if mod(psi,phaseSet{i}) < M.Mrhs.eps
                  if round(biSet{i}((y(1:2))),0) ~= 1
                     psiOut{iPhi} = phaseSet{i};
                     bIOut{iPhi} = biSet{i}; 
                     exitFlag = 1; 
                     break

                  end
               end

           end
           if exitFlag
               break
           end
       end
       phiSet2{iPhi}{1} = t(1:it);
       phiSet2{iPhi}{2} = phi(1:it,:);


    end
else
    phiSet2 = phiSet;
end


end

