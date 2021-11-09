function [theta] = ConvertOToTheta(O)
%CONVERTTHETATOX 

if size(O,1) == 1
       x = zeros(4,size(O,2));
       x(1,:) = cos(O(1,:));
       x(3,:) = sin(O(1,:));
       theta = ConvertXToTheta(x);
else
    error("Dimension not programmed \n")
end

end

