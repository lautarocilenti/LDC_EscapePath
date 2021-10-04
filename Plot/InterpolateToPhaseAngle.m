function [tq,xq] = InterpolateToPhaseAngle(t,x,M)
%INTERPOLATETOPHASEANGLE 
%     it1 = min(find(mod(t,M.Mrhs.T)<M.Mrhs.psiEps));
    
    t2 = linspace(min(t),min(t)+M.Mrhs.T,10000);
    [~,it2] = min(interp1(t,mod(t,M.Mrhs.T),t2));
    
    tq = 0:M.Mrhs.T:max(t);
    iTsmall = find(tq> min(t));
    if ~isempty(iTsmall)
        tq = tq(iTsmall(1):end);
    end
    
%     tq =  t2(it2):M.Mrhs.T:max(t);
    if length(tq)>1
        for d = 1:size(x,2)
            xq(:,d) = interp1(t,x(:,d),tq);
        end
    else
        xq = x(1,:);
    end
   



end

