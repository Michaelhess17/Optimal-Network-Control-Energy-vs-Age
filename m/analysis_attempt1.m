clc
clear all
close all
%% Load data
load('../output/energy.mat')
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
nodeAvg = zeros(size(energy,1));
for i=1:length(energy)
    x = energy(:,:,:,i); % look at one subject
    xr = reshape(x,[size(x,1),max(ci)*max(ci)]);
    xr = bsxfun(@rdivide,log10(xr),nanmax(tiedrank(xr))); %Divide by number of nodes to standardize
        for j=1:length(xr)
        nodeAvg(j, i) = nanmean(xr(j,:)); % store values for all subjects and nodes
        end
end
nodeAvgAvg = mean(nodeAvg,2);

