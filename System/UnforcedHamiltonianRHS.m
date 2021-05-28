function [rhs] = UnforcedHamiltonianRHS(t,x,M)
rhs = [x(2); ...
        -M.a1*x(1)-M.nu*x(2)-M.a3*x(1).^3+x(4); ...
        x(4)*(3*M.a3*x(1)^2+M.a1); ...
        M.nu*x(4)-x(3)];
end
