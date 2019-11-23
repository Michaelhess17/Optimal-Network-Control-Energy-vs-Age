clc
clear all
close all
%% Load data
load('../output/energy.mat')
load('../mat/age.mat')
load('../mat/nki_mats.mat', 'SC', 'coor', 'ci')
load('../output/nodeAvg.mat')
load('../output/nodeAvgGenerated.mat')
load('../output/scrambled_nets.mat')
%%
keep = ~isnan(age);
E = nodeAvg(:,keep)';
D = squeeze(sum(sum(SC(:,:,keep) > 0)));
EGen = nodeAvgGenerated; 
DGen = squeeze(sum(sum(scrambled_nets)));