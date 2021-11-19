function [rhs] = ThreeDuffingHamiltonianRHS(t,x,M)
%
% x(1) = x1
% x(2) = v1
% x(3) = x2
% x(4) = v2
% x(5) = x3
% x(6) = v3
% x(7) = px1
% x(8) = pv1
% x(9) = px2
% x(10) = pv2
% x(11) = px3
% x(12) = pv3
% 
%
fcos = M.F*cos(M.w*t);
rhs = [ x(2); ...
        -M.a1*x(1)-M.nu*x(2)-M.a3*x(1).^3+fcos-M.kc*(2*x(1)-x(3)-x(5))-x(8); ...
        x(4); ...
        -M.a1*x(3)-M.nu*x(4)-M.a3*x(3).^3+fcos-M.kc*(2*x(3)-x(1)-x(5))-x(10); ...
        x(6); ...
        -M.a1*x(5)-M.nu*x(6)-M.a3*x(5).^3+fcos-M.kc*(2*x(5)-x(3)-x(1))-x(12); ...
        ...
        x(8)*(3*M.a3*x(1)^2+M.a1+2*M.kc)-M.kc*x(10)-M.kc*x(12); ...
        -x(7)+M.nu*x(8); ...
         x(10)*(3*M.a3*x(3)^2+M.a1+2*M.kc)-M.kc*x(8)-M.kc*x(12); ...
         -x(9)+M.nu*x(10); ... 
      x(12)*(3*M.a3*x(5)^2+M.a1+2*M.kc)-M.kc*x(8)-M.kc*x(10); ... 
      -x(11)+M.nu*x(12)];
end

