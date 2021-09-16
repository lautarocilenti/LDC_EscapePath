function [Basins] = LoadBasins(M)

folder = "Data\Basins";
if strcmp(M.rhsString,'Unforced') & M.Mrhs.a1 == -1 & M.Mrhs.a3 == 1 & M.Mrhs.nu == 1
    data = load(fullfile(folder,"UnforcedDuffingBasinOverDamped.mat"));
elseif strcmp(M.rhsString,'Unforced') & M.Mrhs.a1 == -1 & M.Mrhs.a3 == .2 & M.Mrhs.nu == .1
    data = load(fullfile(folder,"UnforcedDuffingBasinLarge.mat"));
elseif strcmp(M.rhsString,'Duffing') & M.Mrhs.a1 == 1 & M.Mrhs.a3 == .3 & M.Mrhs.nu == .1 & M.Mrhs.F == .4 & M.Mrhs.w == 1.4
    data = load(fullfile(folder,"Duffing_BasinNumbers_738258.499042.mat"));
else
     Basins = struct([]);
     return
end


Basins = data.Basins;

end

