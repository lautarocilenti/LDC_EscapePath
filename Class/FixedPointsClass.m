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
        names
        iA
        fA
        nFP_Original
        L2;
    end
    
    methods
        function o = FixedPointsClass(fp,phi,stability,solution,iAString,fAString,L2)
            % Constructor
             fpNorm = vecnorm(fp,2,2);
            [~,is] = sort(fpNorm,'ascend');
            o.FP = fp(is,:);
            o.Phi = phi(is);
            o.Stability = stability(is);
            o.Solution = solution(is);
            o.D = size(fp,2);
            o.nFP = size(fp,1);
            o.nFP_Original = o.nFP;
            o.nSaddles = sum(o.Stability<0);
            o.iSaddles = find(o.Stability<0);
            o.L2max = 1.5*max(vecnorm(fp,2,2));
            o.L2 = L2;
            PrepareInterpolant(o);
            CallClassifyFixedPoints(o);
            SetInitalAttractor(o,iAString);
            SetFinalAttractor(o,fAString);
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
                if iFP > o.nFP_Original
                    fp(iFP,:) = o.FP(iFP,:);
                else
                    for d = 1:o.D
                        fp(iFP,d) = o.interpolant{d,iFP}(t);
                    end
                end
            end
            stability = o.Stability;
        end
        
        function [] = CallClassifyFixedPoints(o)
            stability = o.Stability;
            fp = o.FP;
            
            [names] = ClassifyFixedPoints(stability,fp,o.L2);
            o.names = names;
            
        end
        
        function [iA] = SetInitalAttractor(o,iAString)
            for i = 1:length(o.names)
                if contains(o.names{i},iAString)
                    iA = i;
                    o.iA = iA;
                    return
                end
            end
            error("could not find %s",iAString) 
        end
        
                
        function [fA] = SetFinalAttractor(o,fAString)
            for i = 1:length(o.names)
                if contains(o.names{i},fAString)
                    fA = i;
                    o.fA = fA;
                    return
                end
            end
            fA = 0;
            o.fA = 0;
        end
    end
end

