clc
clear all
close all
%% Load data

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
%% Calculate average node energy for all subjects
nodeAvgAD = zeros(size(energyAD,1));
for i=1:length(energyAD)
    x = energyAD(:,:,:,i); % look at one subject
    xr = bsxfun(@rdivide,tiedrank(x),nanmax(tiedrank(x))); %Divide by number of nodes to standardize
        for j=1:length(xr)
            for k = 1:max(ci)
                for l = 1:max(ci)
                    nodeAvgAD(j,i,k,l) = xr(j,k,l); % store values for all subjects and nodes
                end
            end
        end
        nodeAvgAvgAD = squeeze(mean(nodeAvgAD,2));
end
%%
deg = squeeze(sum(SC > 0, 2));
for i = 1:max(ci)
    for j = 1:max(ci)
        k = j + 9*(i-1);
        subplot(9,9,k)
        axis([0, 150, 0, 1])
        plot(deg, squeeze(nodeAvgAD(:,:,i,j)),'k.')
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
    end
end
%%
saveas(gcf, '../figs/png/degree_figures_all_transitions_AD.png')
