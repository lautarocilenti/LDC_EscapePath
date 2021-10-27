function [thetaNew,sNew,State] = EvaluateNelderMead(thetaNew,sNew,thetaCurrent,sCurrent,State,M);
%EVALUATENELDERMEAD 
n = M.dim+1;
d = M.dim;
L = size(thetaCurrent,2);
Ln = L/n;
if M.xcoordinates | floor(Ln) ~= Ln
    error("Dimension issue in Nelder Mead algorithm")
end

simplexState = State.simplexState;
thetaPrev = State.thetaPrev;
sPrev = State.sPrev;

for i = 1:Ln
    j1 = (i-1)*n+1; j2 = (i)*n;
    if simplexState(j1) > 4
        continue
    end
    
    simplex = thetaCurrent(:,j1:j2);
    simplexS = sCurrent(j1:j2);
    [sSorted,is] = sort(simplexS,"ascend");
    iSimplexNew = is(end);
    sSimplexNew = sNew(iSimplexNew);
    
    flag = 0; % 0 - accept new, 1 - reject, 2 - change to prev step
    if simplexState(j1) == 1
        if sSorted(1) < sSimplexNew & sSimplexNew < sSorted(end-1)
            simplexState(j1:j2) = 1; % call standard reflection
            flag = 0;
        elseif sSimplexNew < sSorted(1)
            simplexState(j1:j2) = 2; % call extension
            flag = 1;
        elseif sSimplexNew >= sSorted(end-1)
            if sSimplexNew > sSorted(end)
                simplexState(j1:j2) = 4; % call inside contraction
                flag = 1;
            else
                simplexState(j1:j2) = 3; % call outside contraction
                flag = 1;
            end
        end
    elseif simplexState(j1) == 2 %evaluate extension
        simplexState(j1:j2) = 1;
        if sSimplexNew > sPrev(iSimplexNew);
            flag = 2;
        else
            flag = 0;
        end
    elseif simplexState(j1) == 3 %evaluate outside contraction
        if sSimplexNew <= sPrev(iSimplexNew);
            simplexState(j1:j2) = 1;
            flag = 0;
        else
            simplexState(j1:j2) = 5; %call shrinking
            flag = 1;
        end
    elseif simplexState(j1) == 4 %evaluate inside contraction
        if sSimplexNew <= sSorted(end)
            simplexState(j1:j2) = 1;
            flag = 0;
        else
            simplexState(j1:j2) = 5; %call shrinking
            flag = 1;
        end
    elseif simplexState(j1) == 5  
        flag =0;
    end
    
    if flag == 2
         thetaNew(:,j1:j2) = thetaPrev(:,j1:j2) ;
         sNew(j1:j2) = sPrev(j1:j2);
    elseif flag == 1
        thetaNew(:,j1:j2) = thetaCurrent(:,j1:j2);
        sNew(j1:j2) = sCurrent(j1:j2);
        thetaPrev(:,j1:j2) = thetaNew(:,j1:j2);
        sPrev = sNew(j1:j2);
    else
        thetaPrev(:,j1:j2) = thetaNew(:,j1:j2);
        sPrev = sNew(j1:j2);
    end
        
end
State.simplexState = simplexState;
State.thetaPrev = thetaPrev;
State.sPrev = sPrev;


end

