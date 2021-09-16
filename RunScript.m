AddAllPaths()
ppList = linspace(0,2*pi,21);
ppList = ppList(1:end-1);
save(sprintf("Data/InitialPhase/Phase%2.2f.mat",ppList(1)),'-v7.3');
for i = 1:length(ppList)
    fprintf("Phase List index: %d \n \n",i)
    pp = ppList(i)
    output = Main({'pp'},{pp});
    save(sprintf("Data/InitialPhase/Phase%.2f.mat",ppList(i)),'-v7.3');
end

