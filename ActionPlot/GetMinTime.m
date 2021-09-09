function [tMin] = GetMinTime(phiSet,S)

[minS,iMinS] = max(S);

tMin = max(phiSet{iMinS}{1});
% tMin = 1.01*tMin;
% 

% [maxS,iMaxS] = max(S);
% phi = phiSet{iMaxS}{2};
% t = phiSet{iMaxS}{1};
% p = vecnorm(phi(:,3:4),2,2);
% iStart = min(find(p>.1));
% [pMin,iNewpMin] = min(p(iStart:end));
% tMin = t(iStart+iNewpMin);
% [pMin,ipMin]
% 
% for i = 1:length(phiSet)
%    phi = phiSet{i}{2};
%    t = phiSet{i}{1};
%    p = vecnorm(phi(:,3:4),2,2);
%    
%    dp = diff(p);
%    pMean = movmean(p,100);
%    dp = diff(pMean);
%     
% end


end