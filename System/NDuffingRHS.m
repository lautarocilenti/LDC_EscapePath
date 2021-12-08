function [rhs] = NDuffingRHS(t,x,M)
   L = length(x);
   ix = 1:2:L;
   ixPrev = [L-1,ix(1:end-1)];
   ixNext = [ix(2:end),1]; 
   iv = 2:2:L;
   rhs = zeros(L,1);
   rhs(ix) = x(iv);
   rhs(iv) = -M.a1.*x(ix)-M.nu.*x(iv)-M.a3.*x(ix).^3+M.F.*(cos(M.w.*t))-M.kc.*(2*x(ix)-x(ixPrev)-x(ixNext));
end

