
clc
clear all
close all
%% Load data
load('../output/energyAD.mat')
load('../output/ciprime.mat')
load('../output/nodeAvgAD.mat')
load('../mat/age.mat')
load('../mat/nki_mats.mat', 'coor', 'ci', 'SC')
%%
strength = repmat(tiedrank(squeeze(sum(SC))),[1,1,9,9]);
for i = 1:size(SC,1)
    for j =1:max(ci)
        for k = 1:max(ci)
            keep = strength(:,:,j,k) == i;
            E = nodeAvgAD(:,:,j,k);
            meanE(i,j,k) = nanmean(E(keep));
            stdE(i,j,k) = nanstd(E(keep));
        end
    end    
end
%%
for i = 1:size(SC,3)
    for j = 1:max(ci)
        for k = 1:max(ci)
            nA(:,i,j,k) = nodeAvgAD(:,i,j,k) - meanE(strength(:,i,j,k));
        end
    end
end
%%
for n = 1:size(SC,1)
    for i = 1:9
        for j = 1:9
            x = squeeze(deg(n,:,i,j));
            x = x(:);
            y = squeeze(nA(n,:,i,j));
            y = y(:);
            keep = ~isnan(y) & ~isnan(x);
            correlation(n,i+9*(j-1)) = corr(x(keep),y(keep));
        end
    end
end
figure; imagesc(correlation)