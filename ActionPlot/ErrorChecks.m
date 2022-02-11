function [] = ErrorChecks(M)
%ERRORCHECKS

if M.methodTest
    if M.includePhase
        error("includePhase should be false for method Test\n")
    end
end
end

