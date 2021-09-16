function [] = PL_MinimumSearch(data)
%PL_MINIMUMSEARCH Summary of this function goes here
%   Detailed explanation goes here
msLog = data.msLog;
legendstr = {};
for i = [length(msLog),1]
%     subplot(length(msLog),1,i)
    theta = msLog{i}{1};
    thetaNorm = vecnorm(theta,2,1);
    [thetaNorm,iSortTheta] = sort(thetaNorm,'ascend');
    S = msLog{i}{2};
    S = S(iSortTheta);
    if i == length(msLog)
        plot(thetaNorm,S,'s:','linewidth',1)
    else
        plot(thetaNorm,S,'x--','linewidth',2)
    end
%     legend(sprintf("iter %d",i-1))
    legendstr{end+1} = sprintf("iter %d",i-1);
    hold on
end
legend(legendstr)
xlabel('$||\theta||$ (radians)','interpreter','latex')
ylabel('$S(q)$','interpreter','latex')
title(sprintf('Global Min: %e',data.minS))

end

