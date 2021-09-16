function [InSaddle] = SolutionInSaddle(t, y, Mrhs)
%SOLUTIONINSADDLE
    y = y';
     for iFP = Mrhs.FixedPoints.iSaddles %saddle points
        [fp,~] = Mrhs.FixedPoints.GetFixedPoint(mod(t,Mrhs.T),iFP);
        if vecnorm(fp-y,2,2)<= Mrhs.rS
            InSaddle = true;
            return; 
        end
     end
    InSaddle = false;
end

