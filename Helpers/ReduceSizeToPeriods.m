function [phiSet] = ReduceSizeToPeriods(phiSet,M)
%DOWNSAMPLEPHISET 
if M.saveMemory == 0
    return
end
for i = 1:length(phiSet)
    phi = phiSet{i};
    

    
    if M.saveMemory == 1;
        i1 = [1,length(phi{1})];
        ifall = [1,length(phi{4})];
        
        phiOut{1} = phi{1}(i1);
        phiOut{2} = phi{2}(i1,:);
        phiOut{3} = phi{3};
        phiOut{4} = phi{4}(ifall);
        phiOut{5} = phi{5}(ifall,:);
        phiOut{6} = phi{6}
        phiSet{i} = phiOut;
    else
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
        phi{6} = phi{6};
        phiSet{i} = phi;
    end
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