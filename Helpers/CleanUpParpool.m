function [] = CleanUpParpool()
%CLEANUPPARPOOL 
if CheckIfCluster()
    poolobj = gcp('nocreate');
    delete(poolobj);
end
end

function [result] = CheckIfCluster()
    if (contains(cd,'lcilenti@umd.edu'))
        result = 1;
    else
        result = 0;
    end
end