function [D] = DistanceFromSaddle(phiSet,M)
%DISTANCEFROMSADDLE 
D = zeros(1,length(phiSet));
for i = 1:length(phiSet)
   phi = phiSet{i};
   t = phi{1};
   x = phi{2};
   tf = t(end);
   xf = x(end,1:M.dim)';
   T = 2*pi/M.Mrhs.w;
   tmod = mod(tf,T);
   iSaddle = find(M.Mrhs.FixedPoints.Stability == -1);
   for j = 1:length(iSaddle)
    Saddles(:,j) = M.Mrhs.FixedPoints.GetFixedPoint(tmod,iSaddle(j));
   end
   d = vecnorm(Saddles - xf,2,1);
   D(i) = min(d);
end

end

