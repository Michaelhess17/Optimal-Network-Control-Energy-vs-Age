%%
clc
clear all
close all
%% Load data
load('../output/energySD.mat')
load('../output/energyAD.mat')
load('../mat/age.mat')
load('../mat/nki_mats.mat', 'coor', 'ci', 'SC')
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
%%
nodeAvgSD = zeros(size(energySD,1));
for i=1:length(energySD)
    x = energySD(:,:,:,i); % look at one subject
    xr = reshape(x,[size(x,1),max(ci)*max(ci)]);
    xr = bsxfun(@rdivide,log10(xr),nanmax(tiedrank(xr))); %Divide by number of nodes to standardize
        for j=1:length(xr)
        nodeAvgSD(j, i) = nanmean(xr(j,:)); % store values for all subjects and nodes
        end
end
nodeAvgAvgSD = nanmean(nodeAvgSD,2);
nodeAvgStdSD = nanstd(nodeAvgSD, [], 2);

for i = 1:663
    degree(i,:) = squeeze(sum(SC(:,:,i) > 0));
    R = rich_club_bu(SC(:,:,i) > 0);
    [~,p] = findpeaks(R);
if ~isempty(p)
    optimalRC(i) = p(length(p));
else
    optimalRC(i) = nan;
end
end
choice = nanmean(optimalRC);
degree = mean(degree);
keep = degree > choice;
figure; boxplot(nodeAvgAvgSD, keep)
figure; boxplot(nodeAvgStdSD, keep)
x = nodeAvgStdSD(keep);
y = nodeAvgStdSD(~keep); 
results(1) = ttest2(x,y); 
x = nodeAvgAvgSD(keep);
y = nodeAvgAvgSD(~keep);
results(2) = ttest2(x,y); 

%%
nodeAvgAD = zeros(size(energyAD,1));
for i=1:length(energyAD)
    x = energyAD(:,:,:,i); % look at one subject
    xr = reshape(x,[size(x,1),max(ci)*max(ci)]);
    xr = bsxfun(@rdivide,log10(xr),nanmax(tiedrank(xr))); %Divide by number of nodes to standardize
        for j=1:length(xr)
        nodeAvgAD(j, i) = nanmean(xr(j,:)); % store values for all subjects and nodes
        end
end
nodeAvgAvgAD = nanmean(nodeAvgAD,2);
nodeAvgStdAD = nanstd(nodeAvgAD, [], 2);
figure; boxplot(nodeAvgAvgAD, keep)
figure; boxplot(nodeAvgStdAD, keep)
x = nodeAvgStdAD(keep);
y = nodeAvgStdAD(~keep); 
results(3) = ttest2(x,y); 
x = nodeAvgAvgAD(keep);
y = nodeAvgAvgAD(~keep);
results(4) = ttest2(x,y); 
%%
save('../output/nodeAvgSD_AD.mat','nodeAvgSD','nodeAvgAvgSD','nodeAvgAD','nodeAvgAvgAD', 'nodeAvgStdAD','nodeAvgStdSD')