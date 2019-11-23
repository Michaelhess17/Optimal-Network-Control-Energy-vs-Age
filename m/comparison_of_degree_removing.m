clc
clear all
close all
%% Load data
load('../output/ciprime.mat')
load('../output/energyAD.mat')
load('../output/ciprime.mat')
load('../output/nodeAvgAD.mat')
load('../mat/age.mat')
load('../mat/nki_mats.mat', 'coor', 'SC')
load('../output/nodeAvgResAD.mat')
%%
for i = 1:max(ci)
    for j = 1:max(ci)
        idx_init = ci == i;
        idx_tgt = ci == j;
        idx_bulk = ~idx_init & ~idx_tgt;
        idx_class(:,i,j) = idx_init+2*idx_tgt+3*idx_bulk;
    end
end
save('../output/idx_class.mat')
%%
deg = squeeze(sum(SC > 0, 2));
deg = repmat(deg, [1,1,9,9]);
for n = 1:3
    for i = 0:size(SC,1)
        for j =1:max(ci)
            for k = 1:max(ci)
                keepNodes = deg(:,:,j,k) == i & idx_class(:,j,k) == n;
                E = nodeAvgAD(:,:,j,k);
                meanE(i+1,j,k,n) = nanmean(E(keepNodes));
                stdE(i+1,j,k,n) = nanstd(E(keepNodes));
            end
        end    
    end
end
%%
for i = 1:size(SC,3)
    for j = 1:max(ci)
        for k = 1:max(ci)
            for n =1:3
                meanEClip = meanE(:,j,k,n);
                keepNodes = idx_class(:,j,k) == n;
                nA(keepNodes,i,j,k) = nodeAvgAD(keepNodes,i,j,k) - meanEClip(deg(keepNodes,i,j,k)+1);
            end
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

for n = 1:size(SC,1)
    for i=1:max(ci)
        for j=1:max(ci)
            x = squeeze(deg(n,:,i,j));
            x = x(:);
            y = squeeze(nodeAvgResAD(n,:,i,j));
            y = y(:);
            keep = ~isnan(y) & ~isnan(x);
            correlation2(n,i+9*(j-1)) = corr(x(keep),y(keep));
        end
    end
end
for n = 1:size(SC,1)
    for i=1:max(ci)
        for j=1:max(ci)
            x = squeeze(deg(n,:,i,j));
            x = x(:);
            y = squeeze(nodeAvgAD(n,:,i,j));
            y = y(:);
            keep = ~isnan(y) & ~isnan(x);
            correlation3(n,i+9*(j-1)) = corr(x(keep),y(keep));
        end
    end
end

for n = 1:size(SC,1)
    for i = 1:9
        for j = 1:9
            x = age;
            x = x(:);
            y = squeeze(nA(n,:,i,j));
            y = y(:);
            keep = ~isnan(y) & ~isnan(x);
            correlation4(n,i+9*(j-1)) = corr(x(keep),y(keep));
        end
    end
end

for n = 1:size(SC,1)
    for i = 1:9
        for j = 1:9
            x = age;
            x = x(:);
            y = squeeze(nodeAvgResAD(n,:,i,j));
            y = y(:);
            keep = ~isnan(y) & ~isnan(x);
            correlation5(n,i+9*(j-1)) = corr(x(keep),y(keep));
        end
    end
end

for n = 1:size(SC,1)
    for i = 1:9
        for j = 1:9
            x = age;
            x = x(:);
            y = squeeze(nodeAvgAD(n,:,i,j));
            y = y(:);
            keep = ~isnan(y) & ~isnan(x);
            correlation6(n,i+9*(j-1)) = corr(x(keep),y(keep));
        end
    end
end

subplot(3,4,1)
imagesc(correlation)

set(gca, 'clim',[-1,1])
title('Corrected by Class / Not Parametric')
subplot(3,4,2)
histogram(correlation(:))
xlim([-1,1])
set(gca,'ytick',[])
set(gca,'ylim',[0 3000])
subplot(3,4,3)
imagesc(correlation4)

set(gca, 'clim',[-1,1])

subplot(3,4,4)
histogram(correlation4(:))
xlim([-1,1])
set(gca,'ytick',[])

subplot(3,4,5)
imagesc(correlation2)
set(gca, 'clim',[-1,1])
title('Not Corrected by Class / Parametric')

subplot(3,4,6)
histogram(correlation2(:))
xlim([-1,1])
set(gca,'ytick',[])
set(gca,'ylim',[0 3000])

subplot(3,4,7)
imagesc(correlation5)
set(gca, 'clim',[-1,1])

subplot(3,4,8)
histogram(correlation5(:))
xlim([-1,1])
set(gca,'ytick',[])
set(gca,'ylim',[0 3000])

subplot(3,4,9)
imagesc(correlation3)
set(gca, 'clim',[-1,1])
title('Uncorrected')

subplot(3,4,10)
histogram(correlation3(:))
xlim([-1,1])
set(gca,'ytick',[])
set(gca,'ylim',[0 3000])

subplot(3,4,11)
imagesc(correlation6)
set(gca, 'clim',[-1,1])

subplot(3,4,12)
histogram(correlation6(:))
xlim([-1,1])
set(gca,'ytick',[])
set(gca,'ylim',[0 3000])
%%