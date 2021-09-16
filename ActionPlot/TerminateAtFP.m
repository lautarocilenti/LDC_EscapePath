function [value, isterminal, direction] = TerminateAtFP(t, y, Mrhs)
y = y';

[fp,stability] = Mrhs.FixedPoints.GetAllFixedPoint(mod(t,Mrhs.T));
distance = vecnorm(fp-repmat(y,size(fp,1),1),2,2);
if stability == 1 %attractor
    value = distance-Mrhs.rA;
else
    value = distance-Mrhs.rS;
end

isterminal = ones(size(value));   % Stop the integration
direction  = zeros(size(value));

end

