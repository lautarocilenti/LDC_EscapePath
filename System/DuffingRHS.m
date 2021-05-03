function [rhs] = DuffingRHS(t,x,M)
%RHS Deterministic Unforced Duffing Oscillator
rhs = [x(2); ...
        -M.a1*x(1)-M.nu*x(2)-M.a3*x(1).^3+M.F*cos(M.w*t)+x(4);
        x(4)*(3*M.a3*x(1)^2+M.a1); ...
        M.nu*x(4)-x(3)];
end

