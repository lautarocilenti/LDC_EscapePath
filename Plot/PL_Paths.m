function [] = PL_Paths(phiSet,M)
%PL_PATHS 


for iPhi = 1:length(phiSet)
    t = phiSet{iPhi}{1};
    phi = phiSet{iPhi}{2};
    tFall = phiSet{iPhi}{4};
    phiFall = phiSet{iPhi}{5};
    if strcmp(M.rhsString,'Unforced');
        plot(phi(:,1),phi(:,2))
    else strcmp(M.rhsString,'Duffing');
        color = [rand() rand() rand()];
        it = find(mod(M.Mrhs.w*t,2*pi)<M.Mrhs.psiEps);
        plot(phi(it,1),phi(it,2),'o-','Color',color);
        if ~strcmp(M.terminateType,"Saddle")
            itF = find(mod(M.Mrhs.w*tFall,2*pi)<M.Mrhs.psiEps);
            plot([phi(it(end),1);phiFall(itF,1)],[phi(it(end),2);phiFall(itF,2)],'x--','Color',color);
        end
    end
    hold on
end




end

