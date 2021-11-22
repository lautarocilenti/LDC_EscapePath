function [names] = ClassifyFixedPoints(stability,fp)
%CLASSIFYFIXEDPOINTS 

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
                    if stability(i) == -1 & fp1 == fp2
                        saddleNorm = fp1;
                    end
                end
                for i = 2:size(fp,1)-1
                    fp1 = vecnorm(fp(i,1:2),2,2);
                    fp2 = vecnorm(fp(i,3:4),2,2);
                    
                    if stability(i) == -1
                        names{i} = sprintf("Saddle %d",c); c = c+1;
                        if fp1>fp2 & fp1>saddleNorm
                            names{i} = "HS Saddle";
                        elseif fp2>fp1 & fp2>saddleNorm
                            names{i} = "SH Saddle";
                        elseif fp1==fp2 
                            names{i} = "SS Saddle";
                        elseif fp1>fp2 & fp2<saddleNorm
                            names{i} = "SL Saddle";
                        elseif fp2>fp1 & fp1<saddleNorm
                            names{i} = "LS Saddle";
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
                    fp1 = vecnorm(fp(i,1:2),2,2);
                    fp2 = vecnorm(fp(i,3:4),2,2);
                    fp3 = vecnorm(fp(i,5:6),2,2);
                    fpNorm(i,:) = [round(fp1,4) round(fp2,4) round(fp3,4)];
                    if stability(i) == -1 & fp1 == fp2  & fp1 == fp3;
                        saddleNorm = fp1;
                    end
                end
                
                fpSaddleNorms = fpNorm(stability==-1,:);
                threshHigh = abs(fpHighNorm - saddleNorm)/2;
                threshLow = abs(fpLowNorm - saddleNorm)/2;
                
                iiHigh = fpNorm> fpHighNorm-threshHigh;
                iiLow = fpNorm < fpLowNorm+threshLow;
                iiSaddle = ~(iiHigh | iiLow);
               
                for i = 1:size(iiSaddle,1)
                    letters(i,iiSaddle(i,:)) = "S";
                    letters(i,iiHigh(i,:)) = "H";
                    letters(i,iiLow(i,:)) = "L";
                    if stability(i) == 1
                        letters(i,4) = " Attractor";
                    else
                        letters(i,4) = " Saddle";
                    end
                end
                namesTemp = strcat(letters(:,1),letters(:,2),letters(:,3),letters(:,4));
                for i = 1:length(namesTemp)
                    names{i} = namesTemp(i);
                end
        

        end
end

