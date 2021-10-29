function [phi] = GenerateFullPath(phi,m)
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
    xq = phi{2};
    xo  = xq(1,:)';
    tspan = phi{1};
    [t,y] =  m.solver(m.HamiltonianRHS, tspan, xo,[opts],m.Mrhs);
    phi{1} = t;
    phi{2} = y;
    
    xq = phi{5};
    xo  = xq(1,:)';
    tspan = phi{4};
    [t,y] =  m.solver(m.RHS, tspan, xo,[opts],m.Mrhs);
    phi{4} = t;
    phi{5} = y;

end