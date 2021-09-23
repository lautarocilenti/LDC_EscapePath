function [] = PL_GridSearch(data)
%PL_DESCENT 

msLog = data.msLog{end};

if length(msLog)<3
    return
end

gs = msLog{3};

if gs.Count == 0
    return
end
sCurrent = gs.sCurrent;
stepSizeNormCell = cellfun(@(x) vecnorm(x,2,1),gs.stepSize,"UniformOutput", false );
stepSizeNorm = zeros(length(stepSizeNormCell),size(sCurrent,2));
for i = 1:length(stepSizeNormCell)
    stepSizeNorm(i,:) = stepSizeNormCell{i};
end
for i = 1:size(sCurrent,2)
   s(:,i) = smallerPriorValue(sCurrent(:,i)); 
end

iteration = [1:size(sCurrent,1)]-1;

subplot(3,1,1)
for i = 1:size(sCurrent,2)
   semilogy(iteration,s(:,i),'x-')
   hold on
end
xlabel('iteration')
ylabel('action')

subplot(3,1,2)
iteration = [1:size(s,1)-1];

for i = 1:size(s,2)
   stairs(iteration,(abs(diff(s(:,i)))),'x-')
   hold on
end

xlabel('iteration')
ylabel('descent improvement')
set(gca, 'YScale', 'log')
subplot(3,1,3)
iteration = [1:size(sCurrent,1)]-1;
for i = 1:size(s,2)
   stairs(iteration,(stepSizeNorm(:,i)),'x-')
   hold on
end

xlabel('iteration')
ylabel('norm descent step size')
set(gca, 'YScale', 'log')
% 
% descentTheta{1} = Descent.thetaCurrent{1};
% for i = 1:length(Descent.thetaNew)
% descentTheta{1+i} =  Descent.thetaNew{i};
% end
% descentS = [Descent.sCurrent(1,:);Descent.sNew];
% descentCount = [1:size(descentTheta,1)]-1;
% descentFD = Descent.fdCost;

% 
% subplot(3,1,1)
% theta = vecnorm(data.theta,2,1);
% 
% [theta, isort] = sort(theta,'ascend');
% plot( theta,data.S(isort),'.:b')
% xlabel('$\theta$ (radians)','interpreter','latex')
% ylabel('$S$','interpreter','latex')
%  hold on
% colors = ['rmgk'];
% for i = 1:length(descentTheta)
%     thetaP = vecnorm
%        ii = find([1 ;diff(descentS(:,i))] ~=0);
%        theta = 
%        plot(descentTheta(ii,i),descentS(ii,i),'o-','linewidth',2,'markersize',2,'color',colors(i))
%        hold on
%     end
% end
% subplot(3,1,2)
% for i = 1:size(descentTheta,2)
%    ii = find([1 ;diff(descentS(:,i))] ~=0);
%    semilogy(descentCount(ii),descentS(ii,i),'o-','linewidth',2,'markersize',2,'color',colors(i))
%    hold on
% end
% xlabel('$iteration','interpreter','latex')
% ylabel('$S$','interpreter','latex')
% 
% subplot(3,1,3)
% for i = 1:size(descentTheta,2)
%    plot(descentCount(2:end),abs(descentFD(:,i)),'o-','linewidth',1,'markersize',2,'color',colors(i))
%    hold on
% end


end

function [y] = smallerPriorValue(x)
    m = Inf;
    y = zeros(size(x));
    for i = 1:length(x)
        if x(i) < m
            m = x(i);
        end
            y(i) = m;
    end

end

