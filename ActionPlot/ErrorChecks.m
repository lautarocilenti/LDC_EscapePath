function [] = ErrorChecks(M)
%ERRORCHECKS

if M.methodTest
    if M.includePhase
        error("ParamError: includePhase should be false for method Test\n")
    end
else
       if (M.dim == 2 & ~strcmp(M.rhsString,"Duffing")) | (M.dim == 4 & ~strcmp(M.rhsString,"TwoDuffing")) | (M.dim == 6 & ~(strcmp(M.rhsString,"ThreeDuffing") | (strcmp(M.rhsString,"NDuffing")))) |  (M.dim > 6 & ~strcmp(M.rhsString,"NDuffing")) 
          error("ParamError: Dimension and rhsstring do not match\n")   
       end
        
        
end

end

