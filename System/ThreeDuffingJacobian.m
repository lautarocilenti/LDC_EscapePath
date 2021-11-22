function [res] = ThreeDuffingJacobian(X,M)
% 

J11 = zeros(3,3);
J12 = eye(3);
J13 = zeros(3,3);
J14 = zeros(3,3);

J21 = @(a1,a3,kc,x1,x3,x5) ...
    [-3*a3*x1^2-1-2*kc, kc, kc; ...
    kc, -2*a3*x3^2, kc; ...
    kc, kc, -3*a3*x5^2];
J22 =  @(delta) -eye(3)*delta;
J23 = zeros(3);
J24 = eye(3);

J31 = @(a3,p2,p4,p6,x1,x3,x5) ...
    [ 6*a3*p2*x1, 0, 0; 0 6*a3*p4*x3 0; 0 0 6*a3*p6*x5];
J32 = zeros(3);
J33 = zeros(3);
J34 = @(a1,a3,kc,x1,x3,x5) -J21(a1,a3,kc,x1,x3,x5);

J41 = zeros(3);
J42 = zeros(3);
J43 = -eye(3);
J44 = @(delta) delta*eye(3);

J = @(a1,a3,kc,delta,x1,x3,x5,p2,p4,p6) ... 
    [J11, J12,J13,J14; ...
    J21(a1,a3,kc,x1,x3,x5), J22(delta), J23,J24; ...
    J31(a3,p2,p4,p6,x1,x3,x5),J32,J33,J34(a1,a3,kc,x1,x3,x5); ...
    J41,J42,J43,J44(delta)];





    
res = J(M.a1,M.a3,M.kc,M.nu,X(1),X(2),X(3),0,0,0); 

end

