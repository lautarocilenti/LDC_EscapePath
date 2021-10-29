function  [thetaNew,State] = DetermineFletcherReevesNextSteps(thetaCurrent,sCurrent,iGradient,fdCost,State,M);
%DETERMINEFLETCHERREEVESSTEPSIZE 
%     State 1 starts it off
%     State 2 does backtracking line search
%     State 3 calculates the new gradient based step

    i1 = State.state == 1 & iGradient;
    i2 = State.state == 2 & iGradient;
    i3 = State.state == 3 & iGradient;
    
    if any(i1)
        State.pk =  ResizePk(-fdCost);
        State.dfk = fdCost;
        State.gamma = ones(1,size(thetaCurrent,2));
        State.counter = zeros(1,size(thetaCurrent,2));
        thetaNew = thetaCurrent(:,iGradient)+State.pk;
        return
    end
    
    if any(i3)
        [State] = GenerateNewPk(i3,fdCost,State)
        State.dfk(:,i3) = fdCost(:,i3);
        State.gamma(i3) = 1;
        State.counter(i3) = 0;
    end
    iNew = i2 | i3;
    thetaNew = thetaCurrent(:,iNew)+State.gamma(iNew).*State.pk(:,iNew);


end

function [State] = GenerateNewPk(iRun,fdCost,State)
    dfk1 = fdCost(:,iRun);
    dfk = State.dfk(:,iRun);
    beta =  TransposeMultiplication(dfk1,dfk1)./TransposeMultiplication(dfk,dfk);
    State.pk(:,iRun) = -dfk1+beta.*State.pk(:,iRun);
    State.pk(:,iRun) = ResizePk(State.pk(:,iRun));
end



function [xTy] = TransposeMultiplication(x,y)
    for i = 1:size(x,2)
       xTy(1,i) = x(:,i)'*y(:,i); 
    end
end

function [pk] = ResizePk(pk)
    pkNorm = vecnorm(pk,2,1); 
    ipk = (pkNorm>1);
    if any(ipk)
        pk(:,ipk) = pk(:,ipk)./pkNorm(ipk);
    end
end