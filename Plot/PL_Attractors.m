function [] = PL_Attractors(attractors)
%PL_ATTRACTORS 

for iA = 1:size(attractors,1)
    plot(attractors(iA,1),attractors(iA,2),'.k','markersize',20)
    hold on
end

end

