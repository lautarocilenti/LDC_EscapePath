function [theta] = ConvertXToTheta(x)
%CONVERTXTOTHETA Convert X coordinates to spherical coordinates.
%Coordinates must be in unit n-sphere
if size(x,1) == 4
       theta = zeros(3,size(x,2));
       d1 = vecnorm(x,2,1);
       theta(1,:) = acos(x(1,:)./d1);
       d2 = vecnorm(x(2:end,:),2,1);
       theta(2,:) = acos(x(2,:)./d2);
       d3 = vecnorm(x(3:end,:),2,1);
       ii = x(4,:)>= 0;
       
       theta(3,ii) = acos(x(3,ii)./d3(ii));
       theta(3,~ii) = 2*pi - acos(x(3,~ii)./d3(~ii));
       
%        theta(1,:) = mod(theta(1,:),pi);
%        theta(2,:) = mod(theta(2,:),pi);
%        theta(3,:) = mod(theta(3,:),2*pi);
elseif size(x,1) == 2
    theta = zeros(1,size(x,2));
    theta(1,:) = atan2(x(2,:),x(1,:)); 
    theta(1,:) = mod(theta(1,:),2*pi);
else
    error("Dimension not programmed\n")
end
end

