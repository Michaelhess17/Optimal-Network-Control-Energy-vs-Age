clc
clear all
close all
%% Load data
load('../output/energyAD.mat')
load('../output/nodeAvgAD.mat')
load('../output/ciprime.mat')
load('../mat/age.mat')
load('../mat/nki_mats.mat', 'coor', 'SC')
%%
polySize = 2;
deg = squeeze(sum(SC > 0, 2));
params = zeros(length(deg),max(ci),max(ci),polySize+1);
nA = zeros(size(SC,1), size(SC,3),max(ci),max(ci));
nAFit = zeros(size(SC,1),size(SC,3),max(ci),max(ci));
nodeAvgResAD = zeros(size(SC,1), size(SC,3),max(ci),max(ci));
%%
for i=1:663
    for j = 1:max(ci)
        for k = 1:max(ci)
            nA(:,i,j,k) = nodeAvgAD(:,i,j,k);
            params(i,j,k,:) = polyfit(deg(:,i), nA(:,i,j,k), polySize);
            nAFit(:,i,j,k) = polyval(squeeze(params(i,j,k,:)), deg(:,i));
            nodeAvgResAD(:,i,j,k) = nA(:,i,j,k) - nAFit(:,i,j,k);
        end
    end
end
%%
save('../output/nodeAvgResAD.mat','nodeAvgResAD')