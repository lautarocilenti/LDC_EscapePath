function [names,stability,fp] = ClassifyFixedPoints(stability,fp,L2)
%CLASSIFYFIXEDPOINTS 
        fpNorms = vecnorm(fp,2,2);
        [~,is] = sort(fpNorms,'ascend');
        fp = fp(is,:);
        stability = stability(is);
        if ~isempty(L2)
            L2 = L2(is); L2 = reshape(L2,length(L2),1);
        end
        if size(fp,2) == 2
               
                names{1} = "Low Attractor";
                names{size(fp,1)} = "High Attractor";
                names{2} = "Saddle";
        elseif size(fp,2) == 4
                names{1} = "LL Attractor";
                names{size(fp,1)} = "HH Attractor";
                c = 1;
                for i = 2:size(fp,1)-1
                    fp1 = vecnorm(fp(i,1:2),2,2);
                    fp2 = vecnorm(fp(i,3:4),2,2);
                    [fp1 fp2]
                    if stability(i) == -1 & round(fp1,3) == round(fp2,3)
                        saddleNorm = fp1;
                    end
                end
                for i = 2:size(fp,1)-1
                    fp1 = vecnorm(fp(i,1:2),2,2);
                    fp2 = vecnorm(fp(i,3:4),2,2);
                    
                    if stability(i) == -1
                        names{i} = sprintf("Saddle %d",c); c = c+1;
                        if fp1>fp2 & fp1>saddleNorm
                            names{i} = "LS Saddle";
                        elseif fp2>fp1 & fp2>saddleNorm
                            names{i} = "SL Saddle";
                        elseif fp1==fp2 
                            names{i} = "SS Saddle";
                        elseif fp1>fp2 & fp2<saddleNorm
                            names{i} = "SH Saddle";
                        elseif fp2>fp1 & fp1<saddleNorm
                            names{i} = "HS Saddle";
                        end
                    else
                        if fp1>fp2
                            names{i} = "HL Attractor";
                        else
                            names{i} = "LH Attractor";
                        end
                    end
                end
        elseif size(fp,2) == 6
                 c = 1;
                fpHighNorm = vecnorm(fp(end,1:2),2,2);
                fpLowNorm = vecnorm(fp(1,1:2),2,2);
                for i = 1:size(fp,1)
                    ifp = 1:2:size(fp,2);
                    for j = 1:length(ifp)
                        fpNorm(i,j) = round(vecnorm(fp(i,ifp(j):ifp(j)+1),2,2),4);
                    end
                    if stability(i) == -1 & all(fpNorm(i,:) == fpNorm(i,1))
                        saddleNorm = fpNorm(i,1);
                        sL2 = L2(i);
                    end
                end
                
                
                A = min(fpNorm,[],2);
                B = fpNorm - A;
                iiHigh = B == max(B,[],2);
                iiLow =  B == min(B,[],2);
                iiSaddle = ~iiHigh & ~iiLow;
                
                jj = repmat((L2>= sL2),1,size(fpNorm,2));
                jjSH = jj & iiLow;
                jjSL = ~jj & iiHigh;
                jjHH = false(size(jj)); jjHH(end,:) = true;
                
                iType = size(iiHigh,2)+1;
                iNorm = iType+1;
                for i = 1:size(iiSaddle,1)
                    letters(i,iiSaddle(i,:)) = "S";
                    letters(i,iiHigh(i,:)) = "H";
                    letters(i,iiLow(i,:)) = "L";
%                     letters(i,iiHighSaddle(i,:)) = "H";
%                     letters(i,iiLowSaddle(i,:)) = "L";
%                     letters(i,iiUnknown(i,:)) = "X";
                    if stability(i) == 1
                        letters(i,iType) = " Attractor";
                    else
                        letters(i,iType) = " Saddle";
                        letters(i,jjSH(i,:)) = "S";
                        letters(i,jjSL(i,:)) = "S";
                    end
                    letters(i,jjHH(i,:)) = "H";
                    letters(i,iNorm) = sprintf(" %.2f",round(L2(i),2));
                end
                namesTemp = letters(:,1);
                for i = 2:size(letters,2)
                    namesTemp = strcat(namesTemp,letters(:,i));
                end
                
                for i = 1:length(namesTemp)
                     names{i} = namesTemp(i);
                end
                if length(find(strcmp(names,"SSS Saddle")))>1
                    error("misclassified")
                end
                
        elseif size(fp,2) >= 8
                c = 1;
                fpHighNorm = vecnorm(fp(end,1:2),2,2);
                fpLowNorm = vecnorm(fp(1,1:2),2,2);
                for i = 1:size(fp,1)
                    ifp = 1:2:size(fp,2);
                    for j = 1:length(ifp)
                        fpNorm(i,j) = round(vecnorm(fp(i,ifp(j):ifp(j)+1),2,2),4);
                    end
                    if stability(i) == -1 & all(fpNorm(i,:) == fpNorm(i,1))
                        saddleNorm = fpNorm(i,1);
                        sL2 = L2(i);
                    end
                end
                
                
                A = min(fpNorm,[],2);
                B = fpNorm - A;
                C = sort(B,2,'ascend');
                for i = 1:size(C,1)
                    uC = unique(C(i,:));
                    maxC(i,1) = uC(end);
                    if length(uC)==3
                       midC(i,1)  = uC(2);
                    else
                        midC(i,1) = uC(1);
                    end
                    minC(i,1) = uC(1);
                end
                
                iiMid = B == repmat(midC,1,size(B,2));
                iiHigh = B == repmat(maxC,1,size(B,2));
                iiLow =  B == repmat(minC,1,size(B,2));
                
                %remove mid values
                for i = 1:size(C,1)
                   midValue = max(fpNorm(i,iiMid(i,:))); 
                   lowValue = max(fpNorm(i,iiLow(i,:))); 
                   highValue =  max(fpNorm(i,iiHigh(i,:))); 
                   dLow = abs(midValue-lowValue);
                   dHigh = abs(midValue-highValue);
                   if dLow == 0 & dHigh == 0 %synchronous mode
                      continue
                   elseif dLow<dHigh 
                       iiLow(i,iiMid(i,:)) = true; %mids are closer to low mode
                   else
                       iiHigh(i,iiMid(i,:)) = true; %mids are closer to high mode
                   end
                end

                
                
                jj = repmat((L2>= sL2),1,size(fpNorm,2));
                jjSH = jj & iiLow;
                jjSL = ~jj & iiHigh;
                jjHH = false(size(jj)); jjHH(end,:) = true;
                

                iType = size(iiHigh,2)+1;
                iNorm = iType+1;
                for i = 1:size(iiHigh,1)
                    letters(i,iiMid(i,:)) = "M";
                    letters(i,iiHigh(i,:)) = "H";
                    letters(i,iiLow(i,:)) = "L";
                    if stability(i) == 1
                        letters(i,iType) = " Attractor";
                        
                    else
                        letters(i,iType) = " Saddle";
                        letters(i,iiMid(i,:)) = "s";
                        letters(i,jjSH(i,:)) = "S";
                        letters(i,jjSL(i,:)) = "S";
                    end
                    letters(i,jjHH(i,:)) = "H";
                    letters(i,iNorm) = sprintf(" %.2f",round(L2(i),2));
                end
                namesTemp = letters(:,1);
                for i = 2:size(letters,2)
                    namesTemp = strcat(namesTemp,letters(:,i));
                end
                
                for i = 1:length(namesTemp)
                     names{i} = namesTemp(i);
                end
                if length(find(strcmp(names,"SSS Saddle")))>1
                    error("misclassified")
                end

        end
end

