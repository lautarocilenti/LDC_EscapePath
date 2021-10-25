function [thetaOut,sOut,iOut] = IdentifyInitialSearchValues(theta,S,nLM)

    if size(theta,1) == 1
        [thetaOut,sOut,iOut] = IdentifyLocalMinima(theta,S,nLM);
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


 %Find local minima
    %Neighbor Differentials
    ds1 = [0 diff(sAug)]; 
    ds2 = [diff(sAug) 0];


    iOut = find(ds1<0 & ds2>0); %index of local minima
    sLM = S(iOut); %energy of local minima

    %Limit to max number of local minima to use
    if length(sLM) > nLM
        [~,iSort] = sort(sLM,'ascend');
        iOut = iOut(iSort);
        iOut = iOut(1:nLM);
        sOut = S(iOut);
    end
    thetaOut = thetaAug(iOut);

    
end