function [theta] = ConvertXToTheta(x)
%CONVERTXTOTHETA Convert X coordinates to spherical coordinates.
%Coordinates must be in unit n-sphere
if size(x,1) == 4
       theta = zeros(3,size(x,2));
       theta(1,:) = acos(x(1,:)); 
       theta(2,:) = acos(x(2,:)./sin(theta(1,:)));
       theta(3,:) = acos(x(3,:)./sin(theta(1,:))./sin(theta(2,:)));
       theta(3,:) = theta(3,:).*sign(x(4,:));
else
    error("Dimension not programmed\n")
end
end

