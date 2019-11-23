%%
close all
clear all
clc
%%
load('../output/ciprime.mat')
load('../output/nodeAvgAD.mat')
load('../output/PhenoDfStuff.mat')
load('../output/idx_class.mat')
load('../mat/age.mat')
load('../mat/nki_mats.mat', 'coor', 'SC')
numNodes = 214;
numSys = 9; 
numSubs = 663;
%%
for j = 1:numSys
    for k = 1:numSys
        for n = 1:3
            if j ~= k
                keepSubs = ~isnan(age);
                keepNodes = idx_class(:,j,k) == n;
                E = nodeAvgAD(keepNodes, keepSubs, j, k);
                rho(keepNodes,j,k) = corr(E', age(keepSubs));
            end
        end
    end
end
%%
rho = reshape(rho, [214,81]);
idx_class = reshape(idx_class, [214,81]);
%%
for iClass = 1:3
    for iNode = 1:214
        idx = idx_class(iNode,:) == iClass;
        vals{iClass}(iNode,:) = rho(iNode,idx);
    end
end
%%
a1 = mean(vals{1},2);
a2 = mean(vals{2},2);
a3 = nanmean(vals{3},2);
%%
results = zeros(numSys,10000);
rho1 = rho(:);
for i = 1:numSys
    if sum(a1(ci==i)<0) > sum(a1(ci==i)>0)
        for j = 1:10000
            perm = rho1(randperm(size(rho1, 1),length(a1(ci==i))));
            randAvg = nanmean(perm);
            if randAvg < nanmean(a1(ci==i))
                results(i,j) = 1;
            end
        end
    else
        for j = 1:10000
            perm = rho1(randperm(size(rho1, 1),length(a1(ci==i))));
            randAvg = nanmean(perm);
            if randAvg > nanmean(a1(ci==i))
                results(i,j) = 1;
            end
        end
    end
end
results2 = zeros(numSys,10000);
for i = 1:numSys
    if sum(a2(ci==i)<0) > sum(a2(ci==i)>0)
        for j = 1:10000
            perm = rho1(randperm(size(rho1, 1),length(a2(ci==i))));
            randAvg = nanmean(perm);
            if randAvg < nanmean(a2(ci==i))
                results2(i,j) = 1;
            end
        end
    else
        for j = 1:10000
            perm = rho1(randperm(size(rho1, 1),length(a2(ci==i))));
            randAvg = nanmean(perm);
            if randAvg > nanmean(a2(ci==i))
                results2(i,j) = 1;
            end
        end
    end
end
results3 = zeros(numSys,10000);
for i = 1:numSys
    if sum(a3(ci==i)<0) > sum(a3(ci==i)>0)
        for j = 1:10000
            perm = rho1(randperm(size(rho1, 1),length(a3(ci==i))));
            randAvg = nanmean(perm);
            if randAvg < nanmean(a3(ci==i))
                results3(i,j) = 1;
            end
        end
    else
        for j = 1:10000
            perm = rho1(randperm(size(rho1, 1),length(a3(ci==i))));
            randAvg = nanmean(perm);
            if randAvg > nanmean(a3(ci==i))
                results3(i,j) = 1;
            end
        end
    end
end
%%
