function [thetaNewGridSearch,iGridSearch] = DetermineNewNonGradientPoints(thetaCurrent,iCurrent,stepSize,M)
%DETERMINENEWNONGRADIENTPOINTS 

    if M.xcoordinates
        [thetaNewGridSearch,iGridSearch] = GridSearchOnX(thetaCurrent,iCurrent,stepSize,M);
    else
        [thetaNewGridSearch,iGridSearch] = GridSearch(thetaCurrent,iCurrent,stepSize,M);
        [thetaNewGridSearch] = ModTheta(thetaNewGridSearch)
    end
    

end

function [thetaNew,iGridSearch] = GridSearch(thetaCurrent,iCurrent,stepSize,M)
   
    d = size(thetaCurrent,1);
    n = size(thetaCurrent,2)
    nNewPoints = d*2;
    thetaNew = zeros(d,nNewPoints*n);
    iGridSearch = zeros(1,nNewPoints*n);
   
    for i = 1:n
        j = (i-1)*(nNewPoints)+1;
        k = (i)*(nNewPoints);
        if M.descent.stochasticGridSearch
            v1 = rand(1,d);
        else
            v1 = zeros(1,d); v1(1) = 1;
        end
        v = [v1'/norm(v1) null(v1)];
        step = stepSize(i)*[v -v];
        thetaNew(:,j:k) = thetaCurrent(:,i)+step;
        iGridSearch(j:k) = iCurrent(i);
    end
end

function [xNew,iGridSearch] = GridSearchOnX(xCurrent,iCurrent,stepSize,M)
   
    d = size(xCurrent,1);
    n = size(xCurrent,2)
    nNewPoints = (d-1)*2;
    xNew = zeros(d,nNewPoints*n);
    iGridSearch = zeros(1,nNewPoints*n);
   
    for i = 1:n
        j = (i-1)*(nNewPoints)+1;
        k = (i)*(nNewPoints);
        p = xCurrent(:,i)';
        vr = (2*p)./norm(p); %radial vector
        vt = null(vr); %tangent vectors
        c = repmat(rand(1,size(vt,2)),size(vt,1),1); 
        v1 = sum(c.*vt,2); v1 = v1'./norm(v1); %random tangent vector
        if M.descent.stochasticGridSearch
            v = [v1' null([vr;v1])]; %orthogonal vectors including random tangent vector
        else
            v = [vr' vt]; %orthogonal vectors not stochastic
        end
        step = stepSize(i)*[v -v];
        xNew(:,j:k) = xCurrent(:,i)+step;
        xNew(:,j:k) = xNew(:,j:k)./vecnorm(xNew(:,j:k),2,1);
        iGridSearch(j:k) = iCurrent(i);
    end
end
