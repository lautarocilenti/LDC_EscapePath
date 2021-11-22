function [L] = ThreeDuffingLagrangian(x,M)
%Three Duffing Lagrangian
L = 1/2.*(vecnorm([x(:,8) x(:,10) x(:,12)],2,2)).^2; 
end

