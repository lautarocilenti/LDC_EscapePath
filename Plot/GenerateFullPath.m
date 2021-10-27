function [phi] = generateFullPath(phi,M)

    xq = phi{2};
    xo  = xq(1,:)';

    %Distributed Paths
    phiSetRaw = IntegrateRHS(xo,M);
    % 

    phiSet = PostProcessTrajectories2(phiSetRaw,M);

    phi = phiSet{1};
end