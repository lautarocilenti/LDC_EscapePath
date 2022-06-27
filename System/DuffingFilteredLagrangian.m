function [L] = DuffingFilteredLagrangian(x,M)
%
L = 1/2.*(1./M.Mrhs.wc.*x(:,6)).^2; 
end

