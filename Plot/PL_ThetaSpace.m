function [] = PL_ThetaSpace(data,type)
if nargin == 1
    type = "action";
end

theta = data.theta;
if strcmp(type,"action")
    cost = data.S;
else
    cost = data.C;
end
ii = find(cost> mean(cost) +3*std(cost));
minTheta = mod(theta(:,data.minPhiIndex),2*pi);
cost(ii) = [];
theta(:,ii) = [];
if size(theta,1) <= 1
    return
elseif size(theta,1) == 2
    x = mod(theta(1,:),2*pi);
    y = mod(theta(2,:),2*pi);

    scatter(x,y,40,cost,'filled')    % draw the scatter plot
    ax = gca;

    cb = colorbar;                                     % create and label the colorbar
    cb.Label.String = 'Action';
    xlabel('Theta 1')
    ylabel('Theta 2')
    return
end
x = mod(theta(1,:),pi);
y = mod(theta(2,:),pi);
z = mod(theta(3,:),2*pi);

scatter3(x,y,z,40,cost,'filled')    % draw the scatter plot
ax = gca;
ax.XDir = 'reverse';
view(-31,14)


minTheta(1:2) = mod(minTheta(1:2),pi);
hold on
plot3(minTheta(1),minTheta(2),minTheta(3),'xr','markerSize',40)
cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'Action';
xlabel('Theta 1')
ylabel('Theta 2')
zlabel('Theta 3')

end