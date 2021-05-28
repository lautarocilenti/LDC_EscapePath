function [value, isterminal, direction] = TerminateAtFPOnceAPeriod(t, y, Mrhs)
y = y';
if mod(Mrhs.w*t,2*pi) <= Mrhs.psiEps 
    for iFP = 1:Mrhs.FixedPoints.nFP
        [fp]= Mrhs.FixedPoints.FP(iFP);
        [stability]= Mrhs.FixedPoints.Stability(iFP);
        value(iFP) = vecnorm(fp-y,2,2)-Mrhs.rA;
        if stability == 1 %attractor
            value(iFP) = vecnorm(fp-y,2,2)-Mrhs.rA; %crosses zero from + when distance <= rA
        else %saddle
            value(iFP) = vecnorm(fp-y,2,2)-Mrhs.rS; %crosses zero from + when distance <= rS
        end
    end
else
    value = [1 1 1]; 
end

isterminal = [1 1 1];   % Stop the integration
direction  = [0 0 0];

end

