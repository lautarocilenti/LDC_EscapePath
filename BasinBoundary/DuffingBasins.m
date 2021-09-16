function [bI] = DuffingBasins(M)
%DUFFINGBASINS Summary of this function goes here
%   Detailed explanation goes here
folder = 'Data/Basins/Duffing';
files = dir(fullfile(folder,'*.mat'));
for i = 1:length(files)
    phase(i) = textscan(files(i).name,"Basins_Phi%f.mat");
    data = load(fullfile(folder,files(i).name));
    Basins = data.Basins;
    x = Basins.x(:,1);
    y = Basins.x(:,2);

    B = ones(size(Basins.B)); %initialize all ones

    jj = find(Basins.B~=Basins.iH); %indeces of basin numbers that are not the initial basin
    B(jj) = 3*ones(size(jj)); %all non initial basins are now the number 3

    z = B;
    A{i} = Basins.A;

    L =sqrt(length(x));


    xGrid = reshape(x,L,L);
    yGrid = reshape(y,L,L);
    zGrid = reshape(z,L,L);

    F{i} = griddedInterpolant(xGrid,yGrid,zGrid);
    
    subplot(4,4,i)
    PL_BasinBoundary(F{i},sprintf("$\\phi$ = %f",phase{i}));
    PL_Attractors(A{i});
end

bI = {F,phase,A};



end

