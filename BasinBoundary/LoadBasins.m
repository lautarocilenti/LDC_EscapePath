function [Basins] = LoadBasins(M)

folder = "Data\Basins";
if strcmp(M.rhsString,'Unforced') & M.Mrhs.a1 == -1 & M.Mrhs.a3 == 1 & M.Mrhs.nu == 1
    data = load(fullfile(folder,"UnforcedDuffingBasinOverDamped.mat"));
elseif strcmp(M.rhsString,'Unforced') & M.Mrhs.a1 == -1 & M.Mrhs.a3 == .2 & M.Mrhs.nu == .1
    data = load(fullfile(folder,"UnforcedDuffingBasinLarge.mat"));
end


Basins = data.Basins;

end

