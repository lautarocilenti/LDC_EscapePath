function [] = PL_ThetaSpace(data)
theta = data.theta;
s = data.S;
x = mod(theta(1,:),pi);
y = mod(theta(2,:),pi);
z = mod(theta(3,:),2*pi);

scatter3(x,y,z,40,s,'filled')    % draw the scatter plot
ax = gca;
ax.XDir = 'reverse';
view(-31,14)

minTheta = mod(theta(:,data.minPhiIndex),2*pi);
minTheta(1:2) = mod(minTheta(1:2),pi);
hold on
plot3(minTheta(1),minTheta(2),minTheta(3),'xr','markerSize',40)
cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'Action';
xlabel('Theta 1')
ylabel('Theta 2')
zlabel('Theta 3')

end