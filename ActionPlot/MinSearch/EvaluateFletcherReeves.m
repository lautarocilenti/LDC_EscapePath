function [thetaNew,sNew,State] = EvaluateFletcherReeves(thetaCurrent,thetaNew,sCurrent,sNew,iRun,State,M)
%EVALUATEFLETCHERREEVES 
%     State 1 starts it off
%     State 2 does backtracking line search
%     State 3 calculates the new gradient based step

    fk = sCurrent;
    fgamma = sNew;
    gamma = State.gamma;
    [dfkTpk] = TransposeMultiplication(State.dfk,State.pk);
    
    rho = .8; eta = .5;  loopMax = 10;
    fprintf("Backtracking iteration %d\n", State.counter)
    
    i3 = (fgamma <= fk + eta*gamma.*dfkTpk) | State.counter > loopMax ; % move on
    i2 = ~i3;
    
    %move on
    State.state(i3) = 3;
    
    %continue backtracking
    State.state(i2) = 2;
    State.counter(i2) = State.counter(i2)+1;
    State.gamma(i2) =  gamma(i2).*rho;
    sNew(i2) = sCurrent(i2);
    thetaNew(:,i2) = thetaCurrent(:,i2);
    State.state
end

function [xTy] = TransposeMultiplication(x,y)
    for i = 1:size(x,2)
       xTy(1,i) = x(:,i)'*y(:,i); 
    end
end