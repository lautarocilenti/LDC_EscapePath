function [rhs] = TwoDuffingRHS(t,x,M)
rhs = [x(2); ...
        -M.a1*x(1)-M.nu*x(2)-M.a3*x(1).^3+M.F*cos(M.w*t)-M.kc*(x(1)-x(3)-M.kc*(x(3)-x(1)));
        x(4); ...
        -M.a1*x(3)-M.nu*x(4)-M.a3*x(3).^3+M.F*cos(M.w*t)-M.kc*(x(3)-x(1))-M.kc*(x(3)-x(1))];
end

