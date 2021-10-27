function [] = PL_XSpace(data)
plots3 = true;
M = data.M;
theta = data.theta;
minTheta = theta(:,data.minPhiIndex);
s = data.S;
ii = find(s> mean(s) +2*std(s));
s(ii) = [];
theta(:,ii) = [];


if M.dim == 4
    if ~plots3

        if ~data.M.xcoordinates
            theta = ConvertThetaToX(theta);
        end
        x = theta(1,:)
        y = theta(2,:)
        z = sqrt(theta(3,:).^2+theta(4,:).^2).*sign(theta(4,:));

        scatter3(x,y,z,40,s,'filled')    % draw the scatter plot
        ax = gca;
        ax.XDir = 'reverse';
        view(-31,14)

        minX = theta(:,data.minPhiIndex);
        minZ = sqrt(minX(3).^2+minX(4).^2)*sign(minX(4));
        hold on
        plot3(minX(1),minX(2),minZ,'xr','markerSize',40)
        cb = colorbar;                                     % create and label the colorbar
        cb.Label.String = 'Action';
        xlabel('X 1')
        ylabel('X 2')
        zlabel('sign(X4)*sqrt(X3^2+X4^2)')
    else
 
        if ~data.M.xcoordinates
            theta = ConvertThetaToX(theta);
            minTheta = ConvertThetaToX(minTheta);
        end

        ii = theta(4,:)>=0;
        x = theta(1,ii);
        y = theta(2,ii) ;
        z = theta(3,ii); ; %sqrt(theta(3,:).^2+theta(4,:).^2).*sign(theta(4,:));
        x2 = theta(1,~ii) ;
        y2 = theta(2,~ii) ;
        z2 = theta(3,~ii);
        

        
        
        

        subplot(2,1,1)

        
        scatter3(x,y,z,40,s(ii),'filled')  
        xlabel('X 1')
        ylabel('X 2')
        zlabel('X 3')
        axis 'equal'
        cb = colorbar;                                     % create and label the colorbar
        cb.Label.String = 'Action';
        title('s3+')
        hold on
        if minTheta(4) >=0
            plot3(minTheta(1),minTheta(2),minTheta(3),'xr','markerSize',40)
        end
        subplot(2,1,2)
        scatter3(x2,y2,z2,40,s(~ii),'filled')% draw the scatter plot
        hold on
        title('s3-')
        axis 'equal'
        if minTheta(4) <0
            plot3(minTheta(1),minTheta(2),minTheta(3),'xr','markerSize',40)
        end


        % minX = theta(:,data.minPhiIndex);
        % minZ = sqrt(minX(3).^2+minX(4).^2)*sign(minX(4));
        % hold on
        % plot3(minX(1),minX(2),minZ,'xr','markerSize',40)
        cb = colorbar;                                     % create and label the colorbar
        cb.Label.String = 'Action';
        xlabel('X 1')
        ylabel('X 2')
        zlabel('X 3')

    end
elseif M.dim == 2
        if ~data.M.xcoordinates
            theta = ConvertThetaToX(theta);
            minTheta = ConvertThetaToX(minTheta);
        end

        x = theta(1,:)
        y = theta(2,:)
        scatter(x,y,40,s,'filled')
        cb = colorbar;
        hold on 
         plot(minTheta(1),minTheta(2),'xr','markerSize',40)
end
end