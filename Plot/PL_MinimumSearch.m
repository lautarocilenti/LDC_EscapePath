function [] = PL_MinimumSearch(data)
%PL_MINIMUMSEARCH Summary of this function goes here
%   Detailed explanation goes here
msLog = data.msLog;
legendstr = {};
M = data.M;
di = ceil(length(msLog)/5);
 counter = 5;
for i = length(msLog):-di:1
   
    subplot(5,1,counter)
    counter = counter-1;
    descent = msLog{i}{3};

%     subplot(length(msLog),1,i)
    theta = msLog{i}{1};
    if M.xcoordinates
        theta = ConvertXToTheta(theta);
    end
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
    if i >1
        thetaCurrent = descent.thetaCurrent{i-1};
        if M.xcoordinates
            thetaCurrent = ConvertXToTheta(thetaCurrent);
        end
         thetaCurrent = vecnorm(thetaCurrent,2,1);
         sCurrent = descent.sCurrent(i-1,:);
         plot(thetaCurrent,sCurrent,'or')
    end
end
% legend(legendstr)
xlabel('$||\theta||$ (radians)','interpreter','latex')
ylabel('$S(q)$','interpreter','latex')
title(sprintf('Global Min: %e',data.minS))

end

