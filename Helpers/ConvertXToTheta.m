function [theta] = ConvertXToTheta(x)
%CONVERTXTOTHETA Convert X coordinates to spherical coordinates.
%Coordinates must be in unit n-sphere
if size(x,1) == 4
       theta = zeros(3,size(x,2));
       theta(1,:) = acos(x(1,:)); 
       theta(2,:) = acos(x(2,:)./sin(theta(1,:)));
       theta(3,:) = acos(x(3,:)./sin(theta(1,:))./sin(theta(2,:)));
       theta(3,:) = theta(3,:).*sign(x(4,:));
       
       theta(1,:) = mod(theta(1,:),pi);
       theta(2,:) = mod(theta(2,:),pi);
       theta(3,:) = mod(theta(3,:),2*pi);
elseif size(x,1) == 2
    theta = zeros(1,size(x,2));
    theta(1,:) = atan2(x(2,:),x(1,:)); 
    theta(1,:) = mod(theta(1,:),2*pi);
else
    error("Dimension not programmed\n")
end
end

