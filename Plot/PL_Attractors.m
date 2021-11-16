function [] = PL_Attractors(attractors,names,offset)
%PL_ATTRACTORS 
if nargin == 2
    offset = 0;
end
if isempty(names)
    names = {'1','2','3','4','5','6','7','8','9'};
end
L = size(attractors,1);
Markers = {'s','x','v','+','*','d','^','>','<'};
SaddleColors = repmat([1 0 0], L,1);
AttractorColor = repmat([0 0 1], L,1);
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

