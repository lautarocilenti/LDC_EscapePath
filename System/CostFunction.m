function [SNew,phiSetNew] = CostFunction(theta,M)
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
        SNew = ((x1 - x2).^2 + (x1).^2 + (x2).^4./10);
        for i = 1:size(theta,2)
           phi{1} = 0;
           phi{2} = theta(:,i);
           phi{3} = 0;
           phi{4} = 1;
           phi{5} = theta(:,i);
           phiSetNew{i} = phi;
        end
    else
         phiSetNew = PostProcessTrajectories2(IntegrateRHS(GenerateInitialConditionsFloquet(theta,M),M),M);
         SNew = IntegrateLagrangian(phiSetNew,M); 
         if M.saveMemory
             [phiSetNew] = ReduceSizeToPeriods(phiSetNew,M);
         end
    end
    
end


