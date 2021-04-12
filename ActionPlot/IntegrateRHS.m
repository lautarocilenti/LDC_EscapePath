function [phiSet] = IntegrateRHS(xoSet,M)
%IntegrateRHS - Uses ode45 to integrate Hamiltonian  from initial conditions until
%maximum time or until TerminateAtBoundaryEvent has a value of zero. 

opts = odeset('RelTol',1e-8,'AbsTol',1e-8,'Events', @TerminateAtBoundaryEvent);


nt = length(M.tspan(1):M.dt:M.tspan(2));
phiSet = {};
tSet =  zeros(nt,M.nIC); 
for ixo = 1:size(xoSet,2)
    xo = xoSet(:,ixo);
    [t,y] =  M.solver(M.RHS, M.tspan, xo,[opts],M.Mrhs);
    phiSet{end+1} = {t,y};
    if (max(t) >= max(M.tspan))
        error("Integration terminated prior to reaching boundary\n")
    end
end

end


function [value, isterminal, direction] = TerminateAtBoundaryEvent(T, y, Mrhs)
value = round(Mrhs.bI(y(1:2)'),0) - 2; %zero when at boundary or origin, one otherwise
isterminal = 1;   % Stop the integration
direction  = 0;
end



