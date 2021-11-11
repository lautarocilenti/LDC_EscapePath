function [status,tout,yout] = IntegrateToFixedPoint(t,xo,Mrhs)
%INTEGRATETOFIXEDPOINT This function is called once a period, It follows
%the trajectory in the original system to an attractor or detects if it is
%at a saddle point.
% The output is a status that contains the following possible values
% 1 - "InitialBasin" - solution enters a sphere within rA of initial attractor
% 2 - "AlternateBasin" - solution enters a sphere within rA of another attractor
% 3 - "Saddle" - solution enters a sphere within rS of a saddle

    if SolutionInSaddle(t, xo, Mrhs)
        status = 2;
        tout = [t];
        y(1,:) = xo;
    else
%         opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events', @TerminateAtFP);
%         tspan = [t,Mrhs.tFall+t];
%         [tout,y,~,~,ie] =  Mrhs.solver(Mrhs.RHS, tspan, xo,[opts],Mrhs);
%         if tout(end) >= Mrhs.tFall+t
%             ie = 0;
%         end
%         status = ie;
        opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
        tout = []; yout = [];
        tsol = t;
        for i = 1:50
            tspan = [tsol(end),Mrhs.tFall+tsol(end)];
           
            [tsol,y] =  Mrhs.solver(Mrhs.RHS, tspan, xo,[opts],Mrhs);
            value = TerminateAtFPStandalone(t(end), y(end,:), Mrhs);
            tout = [tout;tsol(1:end-1)];
            yout = [yout;y(1:end-1,:)];
             xo = y(end,:);
            
            if any(value)
                break
            end
        end
        
%         [tout,y,~,~,ie] =  Mrhs.solver(Mrhs.RHS, tspan, xo,[opts],Mrhs);
        if any(value)
            ie = find(value);
        else
            ie = 0;
%             error("ie = 0")
        end
        status = ie;
    end
    
%     figure()
% 
% fp = Mrhs.FixedPoints; 
% ii = find(tout>tout(end)-(2*pi/Mrhs.w));
% t2 = tout(ii);
% y2 = y(ii,:);
% t2= t2-t2(1);
% 
% hold on
% 
% 
% for i = 1:7
%    sol = fp.Solution{i};
%    plot(sol(:,1),vecnorm(sol(:,2:end),2,2))
%    hold on
%    vecnorm(sol(:,2:end),2,1)
% end
% plot(t2,vecnorm(y2(:,1:end),2,2),':','linewidth',3)

end

