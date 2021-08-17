function [phiSet,msLog] = RunMinSearch(phiSet,msLog,M)
%RUNMINSEARCH 
if length(phiSet)<3
    return
end

nLM = M.MS.nLM;
%loop local minima search
for i = 1:M.MS.maxIter 
    fprintf("\nMinSearch Iter %d \n",i)
%     [phiSet,msLog] = MinSearch1DLocalRegion(phiSet,msLog,nLM,M,[1.5,5.5]);
%     [phiSet,msLog] = MinSearch1D_pNorm(phiSet,msLog,nLM,M);
%     S = msLog{end}{2};
%     [tMin] = GetMinTime(phiSet,S);
%     [S2] = IntegrateLagrangianMinTime(phiSet,tMin,M);
%     msLog{end}{2} = S2;
    [phiSet,msLog,TerminateFlag] = MinSearch1DGradientDescent(phiSet,msLog,nLM,M);
    if TerminateFlag
        break
    end
end
%     S = msLog{end}{2};
%     [tMin] = GetMinTime(phiSet,S);
%     [S2] = IntegrateLagrangianMinTime(phiSet,tMin,M);
%     msLog{end}{2} = S2;

end

