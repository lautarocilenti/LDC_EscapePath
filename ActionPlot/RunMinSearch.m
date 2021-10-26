function [phiSet,msLog] = RunMinSearch(phiSet,msLog,M)
    %RUNMINSEARCH 
    nLM = M.MS.nLM;

    TerminateFlag = false;
    %loop local minima search
    for i = 1:M.MS.maxIter 

        fprintf("\nMinSearch Iter %d \n",i)


        [phiSet,msLog,TerminateFlag] = MinSearchAlgorithm(phiSet,msLog,nLM,M);
        if TerminateFlag
                break
        end



    end
    
end


% currentOscillator = M.descent.oscillatorToOptimize;

% % %     if mod(i,25) == 0
% % % %        SaveToRunTimeFile(msLog,phiSet,M,i)
% % %     end
% %     if M.dim == 2
% %         [phiSet,msLog] = MinSearch1D(phiSet,msLog,nLM,M);
% %         TerminateFlag = false;
% %         %     [phiSet,msLog,TerminateFlag] = MinSearch1DGradientDescent(phiSet,msLog,nLM,M);
% %     elseif M.dim == 4
%     if M.descent.optimizePerOscillator
%         [iOrder] =  max(find(i >=M.MS.maxIter./length(M.descent.oscillatorOrder) .*([1:length(M.descent.oscillatorOrder)]-1)));
%         i >=M.MS.maxIter./length(M.descent.oscillatorOrder) .*([1:length(M.descent.oscillatorOrder)]-1)
%         if isempty(iOrder)
%             break
%         end
%         
%         M.descent.oscillatorToOptimize = M.descent.oscillatorOrder(iOrder);
%         
%         if currentOscillator ~= M.descent.oscillatorToOptimize
%             TerminateFlag = false;
%             currentOscillator = M.descent.oscillatorToOptimize;
%             msLog{end}{3}.newStart = true;
%         elseif TerminateFlag
%             continue
%         end
%         
% %         if iOrder == 2
% %             break
% %         end
%         
%         [phiSet,msLog,TerminateFlag] = PerOscillatorGridSearch(phiSet,msLog,nLM,M);
% 
%         
%     elseif M.descent.GradientDescent
%         [phiSet,msLog,TerminateFlag] = MinSearch3DGradientDescent(phiSet,msLog,nLM,M);
%         if TerminateFlag
%             break
%         end
%     else
%         [phiSet,msLog,TerminateFlag] = MinSearchGridSearch(phiSet,msLog,nLM,M);
%         if TerminateFlag
%             break
%         end
%     end
%         
%     end




function [] = SaveToRunTimeFile(msLog,phiSet,M,i)
dateLog = datenum(datetime('now'));
formatOut = 'yyyy_mm_dd_HH_MM_SS_FFF';
fileName = sprintf("Data/ActionPlot/continue%d.mat",i);
save(fileName,'msLog','phiSet','M','dateLog','-v7.3');
end
