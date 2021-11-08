function [] = PlotGenerator(data)
% %PLOTGENERATOR 
% 
% %Plots
% if data.M.clusterRun
%     return
% end

if data.M.methodTest
    
    PL_TestSearch(data);
    
   return 
end


%Generate Gridded interpolant
f1 = figure();
if data.M.Mrhs.w == 1.4 & data.M.dim == 2
    F = GenerateBasinInterpolant(data.M);
    data.basinInterpolant = F; 
    PL_BasinBoundary(data.basinInterpolant,"");
end
PL_MPEP_2DProject(data);

f2 = figure();
PL_MPEP_2DProject(data,"cost");


% 
f5 = figure();
PL_L2MPEPPath(data);
%  ExportPNG(f5,"L2MPEP")
% 
f6 = figure();
PL_MinimumSearch(data);


%  ExportPNG(f6,"MinSearch")
 
 f7 = figure();
 PL_GridSearch(data);
%  
%  f8 = figure();
% PL_L2AllPath(data)

 f9 = figure();
 if data.M.xcoordinates
    PL_XSpace(data)
 else
     
     PL_ThetaSpace(data)
     figure()
     PL_ThetaSpace(data,"cost")
%      PL_XSpace(data)
 end
 
 



end

