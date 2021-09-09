function [L] = TwoDuffingLagrangian(x,M)
%Two Duffing Lagrangian
L = 1/2.*(vecnorm([x(:,6) x(:,8)],2,2)).^2; 
end

