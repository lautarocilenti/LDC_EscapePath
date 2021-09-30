function [theta] = GetInitialTheta(M);
    fprintf("Getting initial theta ...")

    if M.uniformInX
        rng(0);
        nRandom = M.nRVs;
        eta = randn(M.dim,nRandom);
        etaNorm = vecnorm(eta,2,1); 
        etaNorm = repmat(etaNorm,M.dim,1);
        x = eta./etaNorm;
       [x] = ReduceToNumberofInitialPoints(x,M.nIC); %reduce by arclength to neighbors
       if M.xcoordinates
            theta = x;
       else
            [theta] = ConvertXToTheta(x);
       end


    else
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

    fprintf("Done\n")
end 

function [x] = ReduceToNumberofInitialPoints(x,nIC)
    L = .1;
    while(size(x,2)>nIC)
       [x] = RemoveByArcLength(x,L);
       L = L*1.1; 
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