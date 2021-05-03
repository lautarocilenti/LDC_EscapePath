function [phiSet] = IntegrateRHS(xoSet,M)
%IntegrateRHS - Uses ode45 to integrate Hamiltonian  from initial conditions until
%maximum time or until TerminateAtBoundaryEvent has a value of zero. 
if strcmp(M.rhsString,'Unforced')
    opts = odeset('RelTol',1e-8,'AbsTol',1e-8,'Events', @TerminateAtBoundaryEventUnforced);
elseif strcmp(M.rhsString,'Duffing')
    opts = odeset('RelTol',1e-8,'AbsTol',1e-8,'Events', @TerminateAtBoundaryEventDuffing);
end


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


function [value, isterminal, direction] = TerminateAtBoundaryEventUnforced(t, y, Mrhs)
value = round(Mrhs.bI(y(1:2)'),0) - 2; %zero when at boundary or origin, one otherwise
isterminal = 1;   % Stop the integration
direction  = 0;
end

function [value, isterminal, direction] = TerminateAtBoundaryEventDuffing(t, y, Mrhs)

if sqrt((y(1)-Mrhs.A(1))^2+(y(2)-Mrhs.A(2))^2)<Mrhs.r %in sphere
    value = 1;
else
    if mod(Mrhs.w*t,2*pi) <= Mrhs.eps %once a period
        value = round(Mrhs.bI((y(1:2))'),0) == 1;
    else
        value = 1;
    end
end

isterminal = 1;   % Stop the integration
direction  = 0;
end




