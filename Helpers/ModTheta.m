function [theta] = ModTheta(theta,M)
    if M.methodTest
        return
    end
    if M.includePhase
        phase = mod(theta(end,:),2*pi);
        theta(end,:) = [];
    end
    
    if size(theta,1) == 1
        theta = mod(theta,2*pi);
    elseif size(theta,1) == 3
        theta(3,:) = mod(theta(3,:),2*pi);
        theta(1:2,:) = mod(theta(1:2,:),pi);
    elseif size(theta,1) == 5
        theta(5,:) = mod(theta(5,:),2*pi);
        theta(1:4,:) = mod(theta(1:4,:),pi);
    elseif size(theta,1) == 7
        theta(7,:) = mod(theta(7,:),2*pi);
        theta(1:6,:) = mod(theta(1:6,:),pi);
    elseif size(theta,1) == 9
        theta(9,:) = mod(theta(9,:),2*pi);
        theta(1:8,:) = mod(theta(1:8,:),pi);
    else
        error("mod theta Dim not programmed")
    end
    
    if M.includePhase
        theta = [theta;phase];
    end
end

