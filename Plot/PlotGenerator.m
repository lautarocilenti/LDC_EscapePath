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
PL_MPEP_2DProject(data)





% f4 = figure();
% PL_L2Path(data)
% %  ExportPNG(f4,"L2Paths")
%     F = GenerateBasinInterpolant(data.M);
%     data.basinInterpolant = F; 
% %  PL_BasinBoundary(data.basinInterpolant,"");
% f4 = figure()
% % for i = 1:1
% % PL_Path_2DProject(data,i)
% % hold on
% % end
% 
% f4 = figure()

% 
f5 = figure();
PL_L2MPEPPath(data)
%  ExportPNG(f5,"L2MPEP")
% 
f6 = figure();
PL_MinimumSearch(data);
%  ExportPNG(f6,"MinSearch")
 
 f7 = figure()
 PL_GridSearch(data);
%  
%  f8 = figure();
% PL_L2AllPath(data)

 f9 = figure();
 if data.M.xcoordinates
    PL_XSpace(data)
 else
     PL_ThetaSpace(data)
%      PL_XSpace(data)
 end
 
 



end

