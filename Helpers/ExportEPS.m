function [] = ExportEPS(fig,filename)
    Folder = cd;
    f1.PaperUnits = 'inches';
    f1.PaperPosition = [0 0 6 3];
    PlotsFolder = fullfile(Folder, sprintf('Data/Images/%s.eps',filename));
    print(fig,'-depsc',PlotsFolder)

end

