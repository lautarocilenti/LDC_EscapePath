function [L] = UnforcedLagrangian(x,M)
%UNFORCEDLAGRANGIAN 
% plug in the rhs of the hamiltonian systems q2' as x'', add the velocity
% and position terms q2 and q1 and the resultant lagrangian only depends on
% p2
L = 1/2.*(x(:,4)).^2; 
end

