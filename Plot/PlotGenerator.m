function [] = PlotGenerator(data)
% %PLOTGENERATOR 
% 
% %Plots
if data.M.dim == 2
    %Generate Gridded interpolant
    F = GenerateBasinInterpolant(data.M);
    data.basinInterpolant = F; 
    f1 = figure();
    PL_BasinBoundary(data.basinInterpolant,"");
    PL_Attractors(data.attractors);
    PL_Paths(data.phiSet,data.M)
    ExportPNG(f1,"AllPaths")

    f2 = figure();
    PL_PathEnergy(data)
    ExportPNG(f2,"PathEnergy")
    % 
    % % 
    f3 = figure();
    PL_BasinBoundary(data.basinInterpolant,"");
    PL_Attractors(data.attractors);
    PL_MPEP(data);
    ExportPNG(f3,"MPEP")
end



f4 = figure();
PL_L2Path(data)
ExportPNG(f4,"L2Paths")

f5 = figure();
PL_L2MPEPPath(data)
ExportPNG(f5,"L2MPEP")

f6 = figure();
PL_MinimumSearch(data);


end

