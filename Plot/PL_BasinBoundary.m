function [] = PL_BasinBoundary(basinInterpolant)
%PL_BASIN 

x = linspace(-4,4,500);
y = linspace(-3,3,500);
[xGrid,yGrid] = ndgrid(x,y);

zGrid  = basinInterpolant(xGrid,yGrid);

ii = find(round(zGrid,0)==2);
plot(xGrid(ii),yGrid(ii),'.');
hold on

end

