function [theta] = GetInitialTheta(M);

if M.dim == 2
    dtheta =(2*pi)/(M.nIC);
    theta = [dtheta:dtheta:2*pi]; 
elseif M.dim == 4
    dtheta1 =(2*pi)/(M.nIC);
    dtheta23 =(pi)/(M.nIC);
    theta1GridV = [dtheta1:dtheta1:2*pi];
    theta23GridV = [dtheta23:dtheta23:pi];
    [a,b,c] = ndgrid(theta1GridV,theta23GridV,theta23GridV);
    a = reshape(a,1,size(a,1)*size(a,2)*size(a,3));
    b = reshape(b,1,size(b,1)*size(b,2)*size(b,3));
    c = reshape(c,1,size(c,1)*size(c,2)*size(c,3));
    theta = [a;b;c];
end

end 