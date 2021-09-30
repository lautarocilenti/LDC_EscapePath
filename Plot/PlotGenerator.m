function [] = PlotGenerator(data)
% %PLOTGENERATOR 
% 
% %Plots
if data.M.clusterRun
    return
end

if data.M.dim == 2
    %Generate Gridded interpolant
%     F = GenerateBasinInterpolant(data.M);
%     data.basinInterpolant = F; 
%     
%     f1 = figure();
% %     PL_BasinBoundary(data.basinInterpolant,"");
%     PL_Attractors(data.attractors);
%     PL_Paths(data.phiSet,data.M)
%     axis([-4 4 -4 4])
% %      ExportPNG(f1,"AllPaths")
% 
%     f2 = figure();
%     PL_PathEnergy(data)
%      ExportPNG(f2,"PathEnergy")
    % 
    % % 
    f3 = figure();
%     PL_BasinBoundary(data.basinInterpolant,"");
    PL_Attractors(data.attractors);
    PL_MPEP(data);
    axis([-4 4 -4 4])
%      ExportPNG(f3,"MPEP")
end



% f4 = figure();
% PL_L2Path(data)
% %  ExportPNG(f4,"L2Paths")
%     F = GenerateBasinInterpolant(data.M);
%     data.basinInterpolant = F; 
%  PL_BasinBoundary(data.basinInterpolant,"");
f4 = figure()
for i = 1:1
PL_Path_2DProject(data,i)
hold on
end

f4 = figure()
PL_MPEP_2DProject(data)
% 
f5 = figure();
PL_L2MPEPPath2(data)
%  ExportPNG(f5,"L2MPEP")
% 
f6 = figure();
PL_MinimumSearch(data);
%  ExportPNG(f6,"MinSearch")
 
 f7 = figure()
 PL_GridSearch(data);
%  
 f8 = figure();
PL_L2AllPath(data)

 f9 = figure();
 if data.M.xcoordinates
    PL_XSpace(data)
 else
     PL_ThetaSpace(data)
%      PL_XSpace(data)
 end
 
 



end

