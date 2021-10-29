function [] = PL_TestContour(data)
%PL_TESTCONTOUR Summary of this function goes here
x = linspace(-5,5,200);
y = linspace(-5,5,200);
[xx,yy] = meshgrid(x,y);
for i = 1:size(xx,1)
    theta = [xx(i,:);yy(i,:)];
    [zz(i,:),~] = CostFunction(theta,data.M);;
end
contour(xx,yy,zz,[0 .1:10:130],'linewidth',2)
% color1 = [1:-.01:.65]';
% colors = [zeros(size(color1)) color1   color1*.5];
% colormap(flipud(bone));
c = colorbar;
hold on
if data.M.xcoordinates
    p = nsidedpoly(1000, 'Center', [0 0 ], 'Radius', 1);
    plot(p,'FaceColor', 'w')
end

end


