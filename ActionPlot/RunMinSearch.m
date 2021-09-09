function [phiSet,msLog] = RunMinSearch(phiSet,msLog,M)
%RUNMINSEARCH 
if length(phiSet)<3
    return
end

nLM = M.MS.nLM;
%loop local minima search
for i = 1:M.MS.maxIter 
    fprintf("\nMinSearch Iter %d \n",i)
    SaveToRunTimeFile(msLog,phiSet,M)
    if M.dim == 2
        [phiSet,msLog] = MinSearch1D(phiSet,msLog,nLM,M);
        TerminateFlag = false;
        %     [phiSet,msLog,TerminateFlag] = MinSearch1DGradientDescent(phiSet,msLog,nLM,M);
    elseif M.dim == 4
        [phiSet,msLog,TerminateFlag] = MinSearch3DGradientDescent(phiSet,msLog,nLM,M);
    end
    if TerminateFlag
        break
    end
end


end


function [] = SaveToRunTimeFile(msLog,phiSet,M)
dateLog = datenum(datetime('now'));
formatOut = 'yyyy_mm_dd_HH_MM_SS_FFF';
fileName = sprintf("Data/Continue/continue.mat");
save(fileName,'msLog','phiSet','M','dateLog','-v7.3');
end
