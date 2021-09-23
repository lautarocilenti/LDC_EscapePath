function [] = PL_XSpace(data)
% theta = data.theta;
% s = data.S;
% x = theta(1,:)
% y = theta(2,:)
% z = sqrt(theta(3,:).^2+theta(4,:).^2).*sign(theta(4,:));
% 
% scatter3(x,y,z,40,s,'filled')    % draw the scatter plot
% ax = gca;
% ax.XDir = 'reverse';
% view(-31,14)
% 
% minX = theta(:,data.minPhiIndex);
% minZ = sqrt(minX(3).^2+minX(4).^2)*sign(minX(4));
% hold on
% plot3(minX(1),minX(2),minZ,'xr','markerSize',40)
% cb = colorbar;                                     % create and label the colorbar
% cb.Label.String = 'Action';
% xlabel('X 1')
% ylabel('X 2')
% zlabel('sign(X4)*sqrt(X3^2+X4^2)')

theta = data.theta;
if ~data.M.xcoordinates
    theta = ConvertThetaToX(theta);
end
s = data.S;
ii = theta(4,:)>=0;
x = theta(1,ii) +1;
y = theta(2,ii) ;
z = theta(3,ii); ; %sqrt(theta(3,:).^2+theta(4,:).^2).*sign(theta(4,:));
x2 = theta(1,~ii) -1;
y2 = theta(2,~ii) ;
z2 = theta(3,~ii);

scatter3(x,y,z,40,s(ii),'filled')  
hold on
scatter3(x2,y2,z2,40,s(~ii),'filled')% draw the scatter plot
ax = gca;
ax.XDir = 'reverse';
view(-31,14)

% minX = theta(:,data.minPhiIndex);
% minZ = sqrt(minX(3).^2+minX(4).^2)*sign(minX(4));
% hold on
% plot3(minX(1),minX(2),minZ,'xr','markerSize',40)
cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'Action';
xlabel('X 1')
ylabel('X 2')
zlabel('X 3')
title("s3+ center: (1,0,0),s3- center: (-1,0,0)")


end