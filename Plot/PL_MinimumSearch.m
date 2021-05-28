function [] = PL_MinimumSearch(data)
%PL_MINIMUMSEARCH Summary of this function goes here
%   Detailed explanation goes here
msLog = data.msLog;
legendstr = {};
for i = [length(msLog),1]
%     subplot(length(msLog),1,i)
    theta = msLog{i}{1};
    S = msLog{i}{2};
    if i == length(msLog)
        plot(theta,S,'s:','linewidth',1)
    else
        plot(theta,S,'x--','linewidth',2)
    end
%     legend(sprintf("iter %d",i-1))
    legendstr{end+1} = sprintf("iter %d",i-1);
    hold on
end
legend(legendstr)

end

