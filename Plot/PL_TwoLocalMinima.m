function [] = PL_IterPP()
close all
    D = load('ActionPlot_Duffing_Duffing_2021_06_23_17_04_50_559.mat');
    data = D.data;
    msLog = data.msLog;
    
    PL_L2MPEPPath3(data)
    figure()
    F = GenerateBasinInterpolant(data.M);
    data.basinInterpolant = F; 
    PL_BasinBoundary(data.basinInterpolant,"");
    PL_Attractors(data.attractors);
    PL_MPEP2(data)
    legend('Basin Boundary','FP','FP','FP','theta = 1.5104','theta = 5.5168','location','southeast')
    
end