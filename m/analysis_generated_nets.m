clc
clear all
close all
%% Load data
load('../output/energyGenerated.mat')
load('../mat/age.mat')
load('../mat/nki_mats.mat', 'coor', 'ci')
%% Redefine systems
ciprime = zeros(size(ci));
ciprime(ci >= 1 & ci <= 3) = 1;
ciprime(ci >= 4 & ci <= 6) = 2;
ciprime(ci >= 7 & ci <= 8) = 3;
ciprime(ci == 9) = 4;
ciprime(ci >= 10 & ci <= 11) = 5;
ciprime(ci >= 12 & ci <= 13) = 6;
ciprime(ci == 14) = 7;
ciprime(ci >= 15 & ci <= 16) = 8;
ciprime(ci == 17) = 9;
ci = ciprime;
%% Calculate average node energy for all subjects
nodeAvgGen = zeros(size(energy,1), size(energy, 4), 100);
for p = 1:100
    for i=1:size(energy, 4)
        x = energy(:,:,:,i,p); % look at one subject
        xr = reshape(x,[size(x,1),max(ci)*max(ci)]);
        xr = bsxfun(@rdivide,tiedrank(xr),nanmax(tiedrank(xr))); %Divide by number of nodes to standardize
            for j=1:length(xr)
            nodeAvgGen(j,i,p) = xr(j,:); % store values for all subjects and nodes
            end
    end
end
nodeAvgAvg = mean(nodeAvgGen,2);

