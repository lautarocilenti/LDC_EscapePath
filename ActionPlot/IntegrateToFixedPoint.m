function [status,tout,y] = IntegrateToFixedPoint(t,xo,Mrhs)
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
        opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events', @TerminateAtFP);
        tspan = [t,Mrhs.tFall+t];
        [tout,y,~,~,ie] =  Mrhs.solver(Mrhs.RHS, tspan, xo,[opts],Mrhs);
        if tout(end) >= Mrhs.tFall+t
            ie = 0;
        end
        status = ie;
    end



end

