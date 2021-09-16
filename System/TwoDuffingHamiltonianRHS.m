function [rhs] = TwoDuffingHamiltonianRHS(t,x,M)
%
% x(1) = x1
% x(2) = v1
% x(3) = x2
% x(4) = v2
% x(5) = px1
% x(6) = pv1
% x(7) = px2
% x(8) = pv2
% 
%
rhs = [ x(2); ...
        -M.a1*x(1)-M.nu*x(2)-M.a3*x(1).^3+M.F*cos(M.w*t)-M.kc*(x(1)-x(3))-x(6);
        x(4); ...
        -M.a1*x(3)-M.nu*x(4)-M.a3*x(3).^3+M.F*cos(M.w*t)-M.kc*(x(3)-x(1))-x(8);
        x(6)*(2*M.a3*x(1)^2+M.a1+M.kc)-M.kc*x(8);
        -x(5)+M.nu*x(6);
         x(8)*(2*M.a3*x(3)^2+M.a1+M.kc)-M.kc*x(6);
         -x(7)+M.nu*x(8)];
end

