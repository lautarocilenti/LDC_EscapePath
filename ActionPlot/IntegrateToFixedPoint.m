function [status,tout,y] = IntegrateToFixedPoint(t,xo,Mrhs)
%INTEGRATETOFIXEDPOINT This function is called once a period, It follows
%the trajectory in the original system to an attractor or detects if it is
%at a saddle point.
% The output is a status that contains the following possible values
% 1 - "InitialBasin" - solution enters a sphere within rA of initial attractor
% 2 - "AlternateBasin" - solution enters a sphere within rA of another attractor
% 3 - "Saddle" - solution enters a sphere within rS of a saddle


    opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events', @TerminateAtFP);

    tspan = [t,Mrhs.tFall+t];
    Mrhs.temp.told = t;
    [tout,y,~,~,ie] =  Mrhs.solver(Mrhs.RHS, tspan, xo,[opts],Mrhs);
    if tout(end) >= Mrhs.tFall+t
        ie = 0;
    end
    status = ie;
%     if ie == 0
%         status = "Timeout";
%     elseif ie == Mrhs.iA
%         status = "InitialBasin";
%     elseif Mrhs.FixedPoints.Stability(ie)
%         status = "AlternateBasin";
%     elseif Mrhs.FixedPoints.Stability(ie) == -1
%         status = "Saddle";
%     end
% 
%     if status ~=1 
%         T = 2*pi/Mrhs.w;
%         tq = min(tout):T:max(tout);
%         % tq = tout;
%         x1 = y(:,1);
%         x2 = y(:,2);
%         x1q = interp1(tout,x1,tq);
%         x2q = interp1(tout,x2,tq);
% 
%         hold on
%         FP = Mrhs.FixedPoints.FP;
%         for i = 1:size(FP,1)
%            plot(FP(i,1),FP(i,2),'.','markersize',30) 
%            if Mrhs.FixedPoints.Stability(i)==1
%             p = nsidedpoly(1000, 'Center', [FP(i,1),FP(i,2)], 'Radius', Mrhs.rA);
%            else
%             p = nsidedpoly(1000, 'Center', [FP(i,1),FP(i,2)], 'Radius', Mrhs.rS);  
%            end
%             plot(p)
%         end
%         plot(x1q,x2q)
%         axis([-4 4 -4 4])
%     end

end

