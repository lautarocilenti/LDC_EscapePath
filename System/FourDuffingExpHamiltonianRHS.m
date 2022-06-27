function [rhs] = FourDuffingExpHamiltonianRHS(t,x,M)
%

   L = length(x)/2;
   ix = 1:2:L;
   ixPrev = [L-1,ix(1:end-1)];
   ixNext = [ix(2:end),1]; 
   ixNext2 = [ix(3:end),ix(1:2)]; 
   iv = 2:2:L;
   ipx = L+1:2:2*L;
   ipv = L+2:2:2*L;
   ipvPrev = [ipv(end),ipv(1:end-1)];
   ipvNext = [ipv(2:end),ipv(1)];
   ipvNext2 = [ipv(3:end),ipv(1:2)];
   
   rhs = zeros(2*L,1);
   rhs(ix) = x(iv);
   rhs(iv) = -M.a1.*x(ix)-M.nu.*x(iv)-M.a3.*x(ix).^3+M.F.*(cos(M.w.*t))-M.kc.*(3.*x(ix)-x(ixPrev)-x(ixNext)-x(ixNext2))-x(ipv);
   
   rhs(ipx) = x(ipv).*(3*M.a3.*x(ix).^2+M.a1+3.*M.kc)-M.kc.*x(ipvPrev)-M.kc.*x(ipvNext)-M.kc.*x(ipvNext2); 
   rhs(ipv) = -x(ipx)+M.nu.*x(ipv); 
   

end

