function [x] = ConvertThetaToX(theta)
%CONVERTTHETATOX Convert theta to x coordinates
if size(theta,1) == 3
       x = zeros(4,size(theta,2));
       x(1,:) = cos(theta(1,:));
       x(2,:) = sin(theta(1,:)).*cos(theta(2,:));
       x(3,:) = sin(theta(1,:)).*sin(theta(2,:)).*cos(theta(3,:));
       x(4,:) = sin(theta(1,:)).*sin(theta(2,:)).*sin(theta(3,:));
elseif size(theta,1) == 1
       x = zeros(2,size(theta,2));
       x(1,:) = cos(theta(1,:));
       x(2,:) = sin(theta(1,:));
else
    
    error("Dimension not programmed \n")
end
end

