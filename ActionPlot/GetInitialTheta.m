function [theta] = GetInitialTheta(M);

if M.dim == 2
    dtheta =(2*pi)/(M.nIC);
    theta = [dtheta:dtheta:2*pi]; 
elseif M.dim == 4
    if M.xcoordinates
        rng(0);
        nRandom = 15000;
        eta = randn(M.dim,nRandom);
        etaNorm = vecnorm(eta,2,1); 
        etaNorm = repmat(etaNorm,M.dim,1);
        x = eta./etaNorm;
       [x] = ReduceToNumberofInitialPoints(x,M.nIC); %reduce by arclength to neighbors
       theta = x;
%        [theta] = ConvertXToTheta(x);
%        y(1,:) = cos(theta(1,:));
%        y(2,:) = sin(theta(1,:)).*cos(theta(2,:));
%        y(3,:) = sin(theta(1,:)).*sin(theta(2,:)).*cos(theta(3,:));
%        y(4,:) = sin(theta(1,:)).*sin(theta(2,:)).*sin(theta(3,:));
%        all(round(x,:) == round(y,4));
%        y2 = x(4,:);

    else
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

end 

function [x] = ReduceToNumberofInitialPoints(x,nIC)
    L = .05;
    while(size(x,2)>nIC)
        size(x,2)
       [x] = RemoveByArcLength(x,L);
       L = L*1.05; 
    end

end

function [x] = RemoveByArcLength(x,L)
    ixo = 1;
    while(ixo<size(x,2))
       xo = x(:,ixo);
       arcLength = ArcLengthFromXo(x,xo);
       ii = find(arcLength<=L); ii = ii(ii~=ixo);
       x(:,ii) = [];
       ixo = ixo + 1;
    end

end

function [arcLength] = ArcLengthFromXo(x,xo);
    d = x-xo;
    r= 1;
    dNorm = vecnorm(d,2,1);
    arcAngle = 2*asin(dNorm./(2*r));
    arcLength = r*arcAngle;

end