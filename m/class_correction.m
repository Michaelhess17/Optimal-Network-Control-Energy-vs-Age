clear all
close all
clc
%% Load data
load('../output/ciprime.mat')
load('../output/energyAD.mat')
load('../output/ciprime.mat')
load('../output/nodeAvgAD.mat')
load('../mat/age.mat')
load('../mat/nki_mats.mat', 'coor', 'SC')
load('../output/nodeAvgResAD.mat')
%%
deg = squeeze(sum(SC>0,2));
strength = squeeze(sum(SC));
strength(strength == 0) = nan;
strength = tiedrank(strength);
energyADlog = log10(energyAD);
[numNodes,~,numSubs] = size(SC);
numSys = max(ci);
%% Correct ties at 0 strength
for i = 1:size(SC,3)
    offset = size(SC,1) - max(strength(:,i));
    strength(:,i) = strength(:,i) + offset;
end
%%
idx_class = nan(numNodes,numSys,numSys);
for i =1:numSys
    for j =1:numSys
        idx_class(ci == i,i,j) = 1;
        idx_class(ci == j,i,j) = 2;
    end
end
idx_class(isnan(idx_class)) = 3;
%% Degree calculation
meanEnergy = nan(numNodes + 1,numSys,numSys,max(idx_class(:)));
expectedEnergy = nan(size(energyADlog));
offset = 2;
for i = 1:numSys
    for j = 1:numSys
        fprintf('%i %i\n',i,j);
        e = zeros(numNodes,numSubs);
        for ideg = 0:numNodes
            mask = deg >= (ideg - offset) & deg <= (ideg + offset);
            current_energyADlog = squeeze(energyADlog(:,i,j,:));
            if i ~= j
                for k = 1:max(idx_class(:))
                    current_idx_class = repmat(idx_class(:,i,j) == k,[1,numSubs]);
                    big_mask = mask & current_idx_class;
                    vals = current_energyADlog(big_mask);
                    meanEnergy(ideg + 1,i,j,k) = mean(vals);
                    e(big_mask) = mean(vals);
                end
            end
        end
        expectedEnergy(:,i,j,:) = permute(e,[1,3,4,2]);
    end
end
%% Strength Calculation
meanEnergy = nan(numNodes + 1,numSys,numSys,max(idx_class(:)));
expectedEnergy = nan(size(energyADlog));
offset = 2;
for i = 1:numSys
    for j = 1:numSys
        fprintf('%i %i\n',i,j);
        e = zeros(numNodes,numSubs);
        for ideg = 1:numNodes
            mask = strength >= (ideg - offset) & strength <= (ideg + offset);
            current_energyADlog = squeeze(energyADlog(:,i,j,:));
            if i ~= j
                for k = 1:max(idx_class(:))
                    current_idx_class = repmat(idx_class(:,i,j) == k,[1,numSubs]);
                    big_mask = mask & current_idx_class;
                    vals = current_energyADlog(big_mask);
                    meanEnergy(ideg,i,j,k) = mean(vals);
                    e(big_mask) = mean(vals);
                end
            end
        end
        expectedEnergy(:,i,j,:) = permute(e,[1,3,4,2]);
    end
end
%%
for i = 1:size(SC,3)
    for j = 1:max(ci)
        for k = 1:max(ci)
            for n =1:3
                %meanEClip = meanE(:,j,k,n);
                keepNodes = idx_class(:,j,k) == n;
                nA(keepNodes,i,j,k) = energyADlog(keepNodes,j,k,i) ./ expectedEnergy(keepNodes,j,k,i);
            end
        end
    end
end
%%
for n = 1:size(SC,1)
    for i = 1:9
        for j = 1:9
            x = deg(n,:);
            x = x(:);
            y = squeeze(nA(n,:,i,j));
            y = y(:);
            keep = ~isnan(y) & ~isnan(x);
            correlation(n,i+9*(j-1)) = corr(x(keep),y(keep));
        end
    end
end