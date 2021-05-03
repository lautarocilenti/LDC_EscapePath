function [] = PL_Paths(phiSet,M)
%PL_PATHS 


for iPhi = 1:length(phiSet)
    t = phiSet{iPhi}{1};
    phi = phiSet{iPhi}{2};
    if strcmp(M.rhsString,'Unforced')
        plot(phi(:,1),phi(:,2))
    else strcmp(M.rhsString,'Duffing')
        it = find(mod(M.Mrhs.w*t,2*pi)<M.Mrhs.psiEps);
        plot(phi(it,1),phi(it,2))
    end
    hold on
end


end

