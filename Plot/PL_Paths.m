function [] = PL_Paths(phiSet,M)
%PL_PATHS 


for iPhi = 1:length(phiSet)
    t = phiSet{iPhi}{1};
    phi = phiSet{iPhi}{2};

    if strcmp(M.rhsString,'Unforced');
        plot(phi(:,1),phi(:,2))
    else strcmp(M.rhsString,'Duffing');
        color = [rand() rand() rand()];
        [tq,xq] = InterpolateToPhaseAngle(t,phi,M);
        plot(xq(:,1),xq(:,2),'o-','Color',color);
        if M.plotFall
            tFall = phiSet{iPhi}{4};
            phiFall = phiSet{iPhi}{5};
            [tqFall,xqFall] = InterpolateToPhaseAngle(tFall,phiFall,M);
            plot([xq(end,1);xqFall(:,1)],[xq(end,2);xqFall(:,2)],'x--','Color',color);      
        end
    end
    hold on
end




end

