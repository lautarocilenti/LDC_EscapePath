function [value, isterminal, direction] = TerminateAtSaddleEvent(t, y, Mrhs)
    L = Mrhs.FixedPoints.nSaddles+1;
    value = ones(1,L);
    isterminal = ones(1,L);  % Stop the integration
    direction  = zeros(1,L);
    
    if mod(Mrhs.w*t,2*pi) <= Mrhs.psiEps  %once a period
        if sqrt((y(1)-Mrhs.qo(1))^2+(y(2)-Mrhs.qo(2))^2)>=Mrhs.rA1 %in initial condition sphere
           x = y(1:Mrhs.dim)';
            for iFP = Mrhs.FixedPoints.iSaddles %saddle points
                [fp,~] = Mrhs.FixedPoints.GetFixedPoint(mod(t,Mrhs.T),iFP);
                value(iFP) = vecnorm(fp-x,2,2)-Mrhs.rS; %crosses zero from + when distance <= rS
            end
            if norm(y,2)>Mrhs.FixedPoints.L2max
                value(L) = -1; %L2 out of bounds
            end
        end
    end



end


