function [] = PL_Attractors(attractors,names,offset)
%PL_ATTRACTORS 
if size(attractors,2) > 6
    return
end
if nargin == 2
    offset = 0;
end
if isempty(names)
    names = {'1','2','3','4','5','6','7','8','9','10'};
end
L = size(attractors,1);
Markers = {'s','x','v','+','*','d','^','>','<','s','x','v','+','*','d','^','>','<','s','x','v','+','*','d','^','>','<','s','x','v','+','*','d','^','>','<'};
if L >10
    SaddleColors = repmat([1 0 0], L,1).*(rand(L,1)*.5+.5);
    AttractorColor = repmat([0 0 1], L,1).*(rand(L,1)*.5+.5);
else
    SaddleColors = repmat([1 0 0], L,1);
    AttractorColor = repmat([0 0 1], L,1);
end

for iA = 1:L
    marker = Markers{iA};
    msize = 7;
    lwidth = 2;
    if contains(names{iA},"Saddle")
        color = SaddleColors(iA,:);
%          msize = 13;
    else
        color = AttractorColor(iA,:);
    end
    plot(attractors(iA,1+offset),attractors(iA,2+offset),marker,'markersize',msize,'linewidth',lwidth,'color',color,'markerfacecolor',color,'displayname',names{iA})
    hold on
end



end

