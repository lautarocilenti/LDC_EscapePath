function [value, isterminal, direction] = TerminateAtFP(t, y, Mrhs)
y = y';

for iFP = 1:Mrhs.FixedPoints.nFP
    [fp,stability] = Mrhs.FixedPoints.GetFixedPoint(mod(t,Mrhs.T),iFP);
    if stability == 1 %attractor
        value(iFP) = vecnorm(fp-y,2,2)-Mrhs.rA; %crosses zero from + when distance <= rA
    else %saddle
        value(iFP) = vecnorm(fp-y,2,2)-Mrhs.rS; %crosses zero from + when distance <= rS
    end
end


isterminal = [1 1 1];   % Stop the integration
direction  = [0 0 0];

end

