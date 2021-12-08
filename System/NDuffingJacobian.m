function [res] = NDuffingJacobian(X,M)
% 
a1 = M.a1; a3 = M.a3; kc = M.kc; delta = M.nu;

N = length(X)/4;
ix = 1:N;
ixPrev = [ix(end),ix(1:end-1)];
ixNext = [ix(2:end),ix(1)];
ipv = [3*N+1:4*N];
J11 = zeros(N,N);
J12 = eye(N);
J13 = zeros(N,N);
J14 = zeros(N,N);

J21 = zeros(N,N);
J31 = zeros(N,N);
for i = 1:N
   J21(i,ix(i)) = -3*a3*X(ix(i))^2-a1-2*kc;
   J21(i,ixPrev(i)) = kc;
   J21(i,ixNext(i)) = kc;
   
   J31(i,i) = 6*a3*X(ipv(i))*X(ix(i));
end

J22 =  -eye(N)*delta;
J23 = zeros(N);
J24 = eye(N);

J32 = zeros(N);
J33 = zeros(N);
J34 =  -J21;

J41 = zeros(N);
J42 = zeros(N);
J43 = -eye(N);
J44 = delta*eye(N);

J = ...
    [J11, J12,J13,J14; ...
    J21, J22 J23,J24; ...
    J31,J32,J33,J34; ...
    J41,J42,J43,J44];





    
res = J; 

end

