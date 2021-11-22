function [theta] = ConvertXToTheta(x)
%CONVERTXTOTHETA Convert X coordinates to spherical coordinates.
%Coordinates must be in unit n-sphere
if size(x,1) == 6
       theta = zeros(size(x,1)-1,size(x,2));
       for i = 1:size(x,1)-2
           d = vecnorm(x(i:end,:),2,1);
           theta(i,:) = acos(x(i,:)./d);
       end
       d = vecnorm(x(end-1:end,:),2,1);
       ii = x(end,:)>= 0;
       theta(end,ii) = acos(x(end-1,ii)./d(ii));
       theta(end,~ii) = 2*pi - acos(x(end-1,~ii)./d(~ii));
elseif size(x,1) == 4
       theta = zeros(3,size(x,2));
       d1 = vecnorm(x,2,1);
       theta(1,:) = acos(x(1,:)./d1);
       d2 = vecnorm(x(2:end,:),2,1);
       theta(2,:) = acos(x(2,:)./d2);
       d3 = vecnorm(x(3:end,:),2,1);
       ii = x(4,:)>= 0;
       
       theta(3,ii) = acos(x(3,ii)./d3(ii));
       theta(3,~ii) = 2*pi - acos(x(3,~ii)./d3(~ii));
       
elseif size(x,1) == 2
    theta = zeros(1,size(x,2));
    theta(1,:) = atan2(x(2,:),x(1,:)); 
    theta(1,:) = mod(theta(1,:),2*pi);
else
    error("Dimension not programmed\n")
end
end

