function [] = PL_BasinBoundary(basinInterpolant,titleString)
%PL_BASIN 

x = linspace(-5,5,500);
y = linspace(-5,5,500);
[xGrid,yGrid] = ndgrid(x,y);

zGrid  = basinInterpolant(xGrid,yGrid);

ii = find(round(zGrid,0)==3);
plot(xGrid(ii),yGrid(ii),'.');
hold on
title(titleString,'interpreter','latex')
xlabel('x')
ylabel('$\dot{x}$','interpreter','latex')
end

