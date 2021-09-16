function [isCluster] = ProgressBar(n,title)
%PROGRESSBAR Summary of this function goes here
%   Detailed explanation goes here
    isCluster = 0;
    if CheckIfCluster();
        fprintf('\n %s Started:\n',title);
        isCluster = 1;
    else
        fprintf('\n %s Progress:\n',title);
        n25 = ceil(n/4);
        n50 = ceil(n/2);
        n75 = ceil(n*3/4);
        bar = repmat('.',1,n);
        bar(n50) = "^";
        bar([n25,n75]) = "x";
        bar(n) = "<";
        fprintf(['\n' bar '\n\n']);
    end
end

function [result] = CheckIfCluster()
    if (contains(cd,'lcilenti@umd.edu'))
        result = 1;
    else
        result = 0;
    end
end