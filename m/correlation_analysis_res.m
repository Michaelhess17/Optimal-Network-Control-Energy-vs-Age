clc
clear all
close all
%% Load data
load('../output/nodeAvgResAD.mat')
load('../output/ciprime.mat')
load('../mat/age.mat')
load('../mat/nki_mats.mat','SC')
%%
correlation = zeros(size(nodeAvgResAD,1),max(ci),max(ci));
keep = ~isnan(age);
for i=1:max(ci)
    for j=1:max(ci)
        nA = nodeAvgResAD(:,:,i,j);
        correlation(:,i,j) = corr(nA(:,keep)', deg(keep)');
    end
end

%%
deg = squeeze(sum(SC > 0, 2));
keep = ~isnan(age);
figure;
for i = 1:max(ci)
    for j = 1:max(ci)
        k = j + 9*(i-1);
        subplot(9,9,k)
        plot(age(keep), squeeze(nodeAvgResAD(:,keep,i,j)),'k.')
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        axis([0, 100, 0, 1])
    end
end
shg
%%
figure;
for i = 1:max(ci)
    for j = 1:max(ci)
            k = j + 9*(i-1);
            subplot(9,9,k)
            plot(deg, squeeze(nodeAvgResAD(:,:,i,j)),'b.')
            set(gca,'YTickLabel',[]);
            set(gca,'XTickLabel',[]);
            axis([0, 150, -1.2, 1.2])
    end
end
shg
%%
saveas(gcf, '../figs/png/degree_vs_Residuals_AD.png')