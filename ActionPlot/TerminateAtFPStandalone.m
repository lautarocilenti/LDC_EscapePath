function [value,distanceNorm] = TerminateAtFPStandalone(t, y, Mrhs)

[fp,stability] = Mrhs.FixedPoints.GetAllFixedPoint(mod(t,Mrhs.T));
distanceNorm = vecnorm(fp-repmat(y,size(fp,1),1),2,2);
ii = stability == 1;
threshold = Mrhs.rS*ones(size(distanceNorm));
threshold(ii) = Mrhs.rA;
value = distanceNorm <= threshold;

% fprintf([repmat(' %.2f,',1,length(distanceNorm)),'\n'],distanceNorm)
% [mn,imn] = min(distanceNorm);
% fprintf('%f,%d\n',mn,imn)
% if imn == 48
%     dummy = 1;
% end

% 
% isterminal = ones(size(value));   % Stop the integration
% direction  = zeros(size(value));

end

