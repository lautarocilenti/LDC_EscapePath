function [] = PL_MPEP_2DProject(data,type,quantity)

if nargin == 1
   type = "action"; 
   quantity = 1;
elseif nargin == 2
    quantity = 1;
end

%PL_MPEP 
if strcmp(type,"action")
    costSet = data.S;
else
    costSet = data.C;
end

[~,iC] = sort(costSet,'ascend');

for i = 1:quantity
    PL_Path_2DProject(data,iC(i));
end

end

