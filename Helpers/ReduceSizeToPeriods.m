function [phiSet] = ReduceSizeToPeriods(phiSet,M)
%DOWNSAMPLEPHISET 
if ~M.saveMemory
    return
end
for i = 1:length(phiSet)
    phi = phiSet{i};
    t = phi{1};
    y = phi{2};
    [tq,q] = InterpAtPeriods(M,t,y);
    phi{1} = tq;
    phi{2} = q;

    tf = phi{4};
    yf = phi{5};
    [tq,q] = InterpAtPeriods(M,tf,yf);
    phi{4} = tq;
    phi{5} = q;
    phiSet{i} = phi;
end


end

function [timeinPeriods,x] = InterpAtPeriods(m,t,y)
    if length(t) <=1 
        return 
    end
    timeinPeriods = 0:2*pi/m.Mrhs.w:max(t);
    iTsmall = find(timeinPeriods< min(t));
    if ~isempty(iTsmall)
       timeinPeriods = timeinPeriods(iTsmall(end):end);
    end
    for i = 1:size(y,2)
        x(:,i) = interp1(t,y(:,i),timeinPeriods);
    end

end