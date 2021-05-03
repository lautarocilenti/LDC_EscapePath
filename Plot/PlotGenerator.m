function [] = PlotGenerator(data)
%PLOTGENERATOR 

%Plots
f1 = figure();
PL_BasinBoundary(data.basinInterpolant,"");
PL_Attractors(data.attractors);
PL_Paths(data.phiSet,data.M)

f2 = figure();
PL_PathEnergy(data.S,data.M)

% 
f3 = figure();
iS = PL_MPEP(data);
PL_BasinBoundary(data.bi{iS},"");
PL_MPEP(data);

% PL_Attractors(data.attractors);


end

