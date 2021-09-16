function [] = PL_IterPP()
    D = load('ActionPlot_Duffing_Duffing_2021_06_23_17_04_50_559.mat');
    data = D.data;
    msLog = data.msLog;
    
    for i = 1:length(msLog)
        theta = msLog{i}{1};
        S = msLog{i}{2};
        thetaAug = [theta theta(1)];
        sAug = [S;S(1)]; %adding the first term to the end to complete the circle
        ds1 = [0 ;diff(sAug)]; 
        ds2 = [diff(sAug);0];
        iLM = find(ds1<0 & ds2>0); %index of local minima
        sLM = S(iLM); %energy of local minima
        %Limit to max number of local minima to use
%         nLM  = 4;
%         if length(sLM) > nLM
%             [~,iSort] = sort(sLM,'ascend');
%             iLM = iLM(iSort);
%             iLM = iLM(1:nLM);
%             sLM = S(iLM);
%         end
        
        tLM = thetaAug(iLM);
        iL = find(tLM<2);
        sLMVector(1,i) = min(sLM(iL));
        iH = find(tLM>2);
        sLMVector(2,i) = min(sLM(iH));
%         [~,iSort] = sort(tLM,'ascend');
%         sLMVector(:,i) = sLM(iSort);
%         sLM = sLM(iSort);
%         tLM(iSort)
        
    end
    figure(1)
    plot(sLMVector(1,:))
    hold on
    plot(sLMVector(2,:))
    hold off
    legend([sprintf('theta = %.8e, ',tLM(2)),sprintf('S(100) = %.8e',sLMVector(1,end))],[sprintf('theta = %.8e, ',tLM(end)),sprintf('S(100) = %.8e',sLMVector(2,end))])
    ylabel('Action')
    xlabel('iteration')
    
    figure(2)
    
    dsLMVector(1,:) = abs(diff(sLMVector(1,:)))+1E-10;
    dsLMVector(2,:) = abs(diff(sLMVector(2,:)))+1E-10;
    
    x = 1:length(dsLMVector);
    ii = find(dsLMVector(1,:)>=10^-12);
%     semilogy(x(ii),dsLMVector(1,ii),'b')
    stairs(x(:),dsLMVector(1,:),'b');
    
    hold on 
    ii = find(dsLMVector(2,:)>=10^-12);
    stairs(x(:),dsLMVector(2,:),'r')
    set(gca,'YScale','log');
    axis([1,100,1E-9,1E-2])
    hold off
    xlabel('min search iteration')
    ylabel('change in local minima minimum')
    legend(sprintf('theta = %.8e',tLM(2)),sprintf('theta = %.8e',tLM(end)))
end