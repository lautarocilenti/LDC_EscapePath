function [phiSet,msLog] = RunMinSearch(phiSet,msLog,M)
%RUNMINSEARCH 
if length(phiSet)<3
    return
end

nLM = M.MS.nLM;
%loop local minima search
for i = 1:M.MS.maxIter 
    fprintf("\nMinSearch Iter %d \n",i)
    [phiSet,msLog] = MinSearch1D(phiSet,msLog,nLM,M);
end

end

