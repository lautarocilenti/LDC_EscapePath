function [thetaOut,sOut,iOut] = IdentifyInitialSearchValues(theta,S,nLM,M)

    if size(theta,1) == 1
        [thetaOut,sOut,iOut] = IdentifyLocalMinima(theta,S,nLM);
    elseif contains(M.searchAlgorithm,"Simplex")
        [thetaOut,sOut,iOut] = IdentifySmallValuesSimplex(theta,S,nLM,M);
    else
        [thetaOut,sOut,iOut] = IdentifySmallValues(theta,S,nLM);
    end

end

function [thetaOut,sOut,iOut] = IdentifySmallValues(theta,S,nLM) 
    [~,is] = sort(S,'ascend');
    if length(is) > nLM
        is = is(1:nLM);
    end
    thetaOut = theta(:,is);
    sOut = S(is);
    iOut = is;

end

function [thetaOut,sOut,iOut] = IdentifyLocalMinima(theta,S,nLM)

    thetaAug = [theta theta(:,1)];
    sAug = [S S(1)];
 %Find local minima
    %Neighbor Differentials
    ds1 = [0 diff(sAug)]; 
    ds2 = [diff(sAug) 0];


    iOut = find(ds1<0 & ds2>0); %index of local minima
    sOut = S(iOut); %energy of local minima

    %Limit to max number of local minima to use
    if length(sOut) > nLM
        [~,iSort] = sort(sOut,'ascend');
        iOut = iOut(iSort);
        iOut = iOut(1:nLM);
        sOut = S(iOut);
    end
    thetaOut = thetaAug(iOut);

    
end

function [thetaOut,sOut,iOut] = IdentifySmallValuesSimplex(theta,S,nLM,M) 
    n = M.dim+1;
    if size(theta,2) <n
        error("Not enough initial points to build simplex\n")
    end
    NLM = floor(size(theta,2)./n);
    if NLM> nLM
        nLM = NLM;
    end
    [~,is] = sort(S,'ascend');
    if length(is) > nLM
        is = is(1:nLM);
    end
    thetaLowS = theta(:,is);
    
    for i = 1:size(thetaLowS,2)
        j1 = (i-1)*n+1; j2 = (i)*n;
        d = vecnorm(theta-thetaLowS(:,i),2,1);
        [~,id] = sort(d,'ascend');
        thetaOut(:,j1:j2) = theta(:,id(1:n));
        sOut(j1:j2) = S(id(1:n));
        iOut(j1:j2) = id(1:n);
    end


end