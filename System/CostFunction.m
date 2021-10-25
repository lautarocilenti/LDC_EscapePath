function [SNew,phiSetNew] = CostFunction(theta,M,saveMemory)
%CostFunction computes the action to escape the basin following a 
%   hamiltonian trajectory with initial position theta. Test
%   mode calculates the cost of a simple 2D function.
%
% Inputs: theta (initial positions), M (struct of parameters)
%         saveMemory (returns trajectories when false)
% Outputs: SNew (action values), phiSetNew (set of trajectories)

    if nargin == 2
        saveMemory = true;
    end
    
    if M.methodTest
        x1 = theta(1,:);
        x2 = theta(2,:);
        SNew = ((x1 - x2).^2 + (x1 - 2).^2 + (x2 - 3).^4). / 10;
        phiSetNew = {};
    else
         phiSetNew = PostProcessTrajectories2(IntegrateRHS(GenerateInitialConditionsFloquet(theta,M),M),M);
         SNew = IntegrateLagrangian(phiSetNew,M); 
         if saveMemory
           phiSetNew = {};
         end
%          [phiSetNew] = ReduceSizeToPeriods(phiSetNew,M);
    end
    
end


