function [contourHandle] = ContourGridVectors(x,y,z,plotType)
%CONTOURGRIDVECTORS Summary of this function goes here
%   Detailed explanation goes here
if nargin == 3
    plotType = "contour";
end


if size(z,1) == size(z,2)
    error("grid vector inputs only in ContourGridVectors\n");
end

L =sqrt(length(x));

if ~(isfinite(L) & L==floor(L))
    error("cannot reshape into square matrix\n");
end


xq = reshape(x,L,L);
yq = reshape(y,L,L);
zq = reshape(z,L,L);


if(plotType == "contour")
    contourHandle = contour(xq,yq,zq);
elseif (plotType == "surf")
    contourHandle = surf(xq,yq,zq);
elseif (plotType == "mesh")
    contourHandle = mesh(xq,yq,zq);
else
    error("unknown plot type, %s",plotType)
end

end

