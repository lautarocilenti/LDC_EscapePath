function [] = PL_Attractors(attractors,offset,names)
%PL_ATTRACTORS 
if nargin == 1
    offset = 0;
end
L = size(attractors,1);
Markers = {'+','*','v','d','^','s','>','<'};
SaddleColors = repmat([1 0 0], L,1);
AttractorColor = repmat([0 0 1], L,1);
for iA = 1:L
    marker = Markers{iA};
    msize = 7;
    if contains(names{iA},"Saddle")
        color = SaddleColors(iA,:);
    else
        color = AttractorColor(iA,:);
    end
    plot(attractors(iA,1+offset),attractors(iA,2+offset),marker,'markersize',msize,'color',color,'markerfacecolor',color,'displayname',names{iA})
    hold on
end



end

