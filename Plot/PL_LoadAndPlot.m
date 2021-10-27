function [] = PL_LoadAndPlot(filepath)
%PL_LOADANDPLOT 

if nargin < 1
    folder = "D:\Melody\OneDrive - University of Maryland\University of Maryland\ResearchData\ActionPlotData\TwoOscillators\October2021";
    filename = "ActionPlot_TwoDuffing_TwoOscillator_ThetaGridSearch_iA4_2021_10_06_15_42_01_243";
    filepath = fullfile(folder,filename);
end

content = load(filepath); 
data = content.data;

fullpath = true;
PlotGenerator(data,fullpath)

end

