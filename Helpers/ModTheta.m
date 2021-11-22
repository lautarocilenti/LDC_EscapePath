function [theta] = ModTheta(theta,M)
    if M.methodTest
        return
    end
    if size(theta,1) == 1
        theta = mod(theta,2*pi);
    elseif size(theta,1) == 2
        theta(1,:) = mod(theta(1,:),2*pi);
        theta(2,:) = mod(theta(2,:),2*pi);
    elseif size(theta,1) == 3
        theta(3,:) = mod(theta(3,:),2*pi);
        theta(1:2,:) = mod(theta(1:2,:),pi);
    elseif size(theta,1) == 5
        theta(5,:) = mod(theta(5,:),2*pi);
        theta(1:4,:) = mod(theta(1:4,:),pi);
    else
        error("mod theta Dim not programmed")
    end
end

