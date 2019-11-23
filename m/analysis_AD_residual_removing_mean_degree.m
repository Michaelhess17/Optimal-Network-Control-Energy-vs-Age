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
deg = repmat(squeeze(sum(SC > 0, 2)),[1,19,9]);
for i = 1:size(SC,1)
    for i =1:9
        for k = 1:9
                keep = deg == i;
                E = nodeAvgAD(:,:,j,k);
                meanE(i,j,k) = nanmean(E(keep));
        end
    end
end
%%
for i = 1:size(SC,3)
   nA(:,i,:,:) = bsxfun(@rdivide,squeeze(nodeAvgAD(:,i,:,:)),meanE');
end
%%;
for i = 1:9
    for j = 1:9
        x = squeeze(deg(:,:,i,j));
        x = x(:);
        y = squeeze(nA(:,:,i,j));
        y = y(:);
        keep = ~isnan(y) & ~isnan(x);
        correlation(i,j) = corr(x(keep),y(keep));
    end
end
imagesc(correlation)
load('../output/nodeAvgResAD.mat')
for i=1:max(ci)
    for j=1:max(ci)
        x = squeeze(deg(:,:,i,j));
        x = x(:);
        y = squeeze(nodeAvgResAD(:,:,i,j));
        y = y(:);
        keep = ~isnan(y) & ~isnan(x);
        correlation2(i,j) = corr(x(keep),y(keep));
    end
end
figure; imagesc(correlation2)