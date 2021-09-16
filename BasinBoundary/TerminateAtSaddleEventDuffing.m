function [value, isterminal, direction] = TerminateAtSaddleEvent(t, y, Mrhs)
    value = [1 1 1];
    isterminal = [1 1 1];   % Stop the integration
    direction  = [0 0 0];
    
    if mod(Mrhs.w*t,2*pi) <= Mrhs.psiEps  %once a period
        if sqrt((y(1)-Mrhs.qo(1))^2+(y(2)-Mrhs.qo(2))^2)>=Mrhs.rA1 %in initial condition sphere
                iiFP = find(1:Mrhs.FixedPoints.Stability<0);
                for iFP = find(1:Mrhs.FixedPoints.Stability<0) %saddle points
                    [fp,~] = Mrhs.FixedPoints.GetFixedPoint(mod(t,Mrhs.T),iFP);
                    value(iFP) = vecnorm(fp-y,2,2)-Mrhs.rS; %crosses zero from + when distance <= rS
                end
                xo = [y(1);y(2)];
                status = IntegrateToFixedPoint(t,xo,Mrhs);
                if status == 0
                    error("Timeout error")
                elseif status >1
                    value(status) = -1;
                end
        end
    end



end


