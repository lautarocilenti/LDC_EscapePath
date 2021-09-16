function [tq,xq] = InterpolateToPhaseAngle(t,x,M)
%INTERPOLATETOPHASEANGLE 
    it1 = min(find(mod(t,M.Mrhs.T)<M.Mrhs.psiEps));
    tq =  t(it1):M.Mrhs.T:max(t);
    if length(tq)>1
        for d = 1:size(x,2)
            xq(:,d) = interp1(t,x(:,d),tq);
        end
    else
        xq = x(1,:);
    end
   



end

