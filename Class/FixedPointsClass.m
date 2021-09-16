classdef FixedPointsClass <handle
    %This handle class is used to refer to large fixed point data without
    %making too many copies of the data
    
    properties
        FP %fixed points at phase angle
        Phi % phase angle
        Stability %stability
        Solution %periodic trajectory
        D %dimension
        nFP %number of fixed points
        L2max
        nSaddles
        iSaddles
        interpolant
    end
    
    methods
        function o = FixedPointsClass(fp,phi,stability,solution)
            % Constructor
            o.FP = fp;
            o.Phi = phi;
            o.Stability = stability;
            o.Solution = solution;
            o.D = size(fp,2);
            o.nFP = size(fp,1);
            o.nSaddles = sum(stability<0);
            o.iSaddles = find(stability<0);
            o.L2max = 1.5*max(vecnorm(fp,2,2));
            PrepareInterpolant(o);
        end
        
        function [] = PrepareInterpolant(o)
            for iFP = 1:o.nFP
                for d = 1:o.D
                   t = o.Solution{iFP}(:,1);
                   tq = linspace(min(t),max(t),2*length(t));
                   xq = interp1(t,o.Solution{iFP}(:,1+d),tq);
                   o.interpolant{d,iFP} = griddedInterpolant(tq,xq);
                end
            end
        end
        
        function [fp,stability] = GetFixedPoint(o,t,iFP)
           %returns interpolated fixed point given time and solution
           fp = zeros(1,o.D);
           for d = 1:o.D
                fp(1,d) = o.interpolant{d,iFP}(t);
           end
           stability = o.Stability(iFP);
        end
        
        function [fp,stability] = GetAllFixedPoint(o,t)
            fp = zeros(o.nFP,o.D);
            for iFP = 1:o.nFP
                for d = 1:o.D
                    fp(iFP,d) = o.interpolant{d,iFP}(t);
                end
            end
            stability = o.Stability;
        end
    end
end

