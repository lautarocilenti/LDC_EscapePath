function [] = PL_Paths(phiSet,M)
%PL_PATHS 

for iPhi = 1:length(phiSet)
    phi = phiSet{iPhi}{2};
    plot(phi(:,1),phi(:,2))
    hold on
end

end

