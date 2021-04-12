function [] = SaveToFile(data,M)
dateLog = datenum(datetime('now'));
formatOut = 'yyyy_mm_dd_HH_MM_SS_FFF';
fileName = sprintf("Data/ActionPlot/ActionPlot_%s_%s_%s.mat",M.rhsString,M.paramNote,datestr(dateLog,formatOut));
save(fileName,'data','dateLog','-v7.3');
end

