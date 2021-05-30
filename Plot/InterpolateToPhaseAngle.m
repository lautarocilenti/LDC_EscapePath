function [tq,xq] = InterpolateToPhaseAngle(t,x,M)
%INTERPOLATETOPHASEANGLE 
    it1 = min(find(mod(t,M.Mrhs.T)<M.Mrhs.psiEps));
    tq =  t(it1):M.Mrhs.T:max(t);
    for d = 1:size(x,2)
        xq(:,d) = interp1(t,x(:,d),tq);
    end



end

