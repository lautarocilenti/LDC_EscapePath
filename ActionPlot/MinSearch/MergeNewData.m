function [thetaOut,sOut,phiSetOut] = MergeNewData(theta,S,phiSet,thetaNewSearch,sNewSearch,phiSetNewSearch)
%MERGENEWDATA
    thetaOut = [theta thetaNewSearch];
    sOut = [S  sNewSearch];
    phiSetOut = [phiSet phiSetNewSearch];
end

