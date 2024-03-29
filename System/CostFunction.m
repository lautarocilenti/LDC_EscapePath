function [costNew,phiSetNew] = CostFunction(theta,M)
%CostFunction computes the action to escape the basin following a 
%   hamiltonian trajectory with initial position theta. Test
%   mode calculates the cost of a simple 2D function.
%
% Inputs: theta (initial positions), M (struct of parameters)
%         
% Outputs: SNew (action values), phiSetNew (set of trajectories)

    
    if M.methodTest
        x1 = theta(1,:);
        x2 = theta(2,:);
        ii = x2>0;
        costNew(~ii) = ((x1(~ii) - x2(~ii)).^2 + (x1(~ii)).^2 + (x2(~ii)).^4./10);
        costNew(ii) = 180*(1/sqrt(2*pi)*exp(1/25.*(-x1(ii).^2-3*x2(ii).^2)))+.1;
        for i = 1:size(theta,2)
           phi{1} = 0;
           phi{2} = theta(:,i);
           phi{3} = 0;
           phi{4} = 1;
           phi{5} = theta(:,i);
           phi{6} = costNew(i);
           phiSetNew{i} = phi;
        end
    else
         phiSetNew = PostProcessTrajectories2(IntegrateRHS(GenerateInitialConditionsFloquet(theta,M),M),M);
         if strcmp(M.costType,"distancefromSaddle")
             [~,phiSetNew] = IntegrateLagrangian(phiSetNew,M);
             [costNew] = DistanceFromSaddle(phiSetNew,M);
         else
             [costNew,phiSetNew] = IntegrateLagrangian(phiSetNew,M);
         end
        
         if M.saveMemory
             [phiSetNew] = ReduceSizeToPeriods(phiSetNew,M);
         end
    end
    
end


