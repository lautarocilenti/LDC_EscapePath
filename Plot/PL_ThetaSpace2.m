function [] = PL_ThetaSpace2(theta,s)

x = mod(theta(1,:),pi);
y = mod(theta(2,:),pi);
z = mod(theta(3,:),2*pi);

scatter3(x,y,z,40,s,'filled')    % draw the scatter plot
ax = gca;
ax.XDir = 'reverse';
% view(-31,14)
% [sSorted,iSortS] = sort(s,'ascend');
% theta2 = [x;y;z];
% thetaSorted = theta2(:,iSortS);

% hold on
% plot3(thetaSorted(1,1:3),thetaSorted(2,1:3),thetaSorted(3,1:3),'xr','markerSize',40)
cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'Action';
xlabel('Theta 1')
ylabel('Theta 2')
zlabel('Theta 3')

end