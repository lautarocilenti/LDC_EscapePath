function [phi] = GenerateFullPath(phi,m)
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);


    xq = phi{2};
    xo  = xq(1,:)';
%     m.xcoordinates = true;
%     m.iA = 7;
%     [~,phiSet] = CostFunction(xo(1:4,:),m);
%     phi = phiSet{1};
%     return
    
    tspan = [phi{1}(1) phi{1}(end)];
    [t,y] =  m.solver(m.HamiltonianRHS, tspan, xo,[opts],m.Mrhs);
    phi{1} = t;
    phi{2} = y;
    
    xq = phi{5};
    xo  = xq(1,:)';
    tspan = [phi{4}(1) phi{4}(end)];
    [t,y] =  m.solver(m.RHS, tspan, xo,[opts],m.Mrhs);
    phi{4} = t;
    phi{5} = y;

end