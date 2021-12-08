function [L] = NDuffingLagrangian(x,M)
%N Duffing Lagrangian
N = size(x,2)/2;
ipv = N+2:2:2*N;
L = 1/2.*(vecnorm([x(:,ipv)],2,2)).^2; 
end

