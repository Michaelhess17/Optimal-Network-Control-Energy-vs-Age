load('C:\Users\Micha\Desktop\Computer Science\MATLAB\nki+controllability\output\nodeAvgGenerated.mat')
load('C:\Users\Micha\Desktop\Computer Science\MATLAB\nki+controllability\output\scrambled_nets.mat')
load('C:\Users\Micha\Desktop\Computer Science\MATLAB\nki+controllability\mat\nki_mats.mat', 'coor', 'ci')
nodeAvgAvg = squeeze(mean(mean(nodeAvgGenerated, 2),3));
nodeAvgStd = squeeze(mean(std(nodeAvgGenerated, [], 2),3));
figure;scatter3(coor(:,1),coor(:,2),coor(:,3), 100, nodeAvgAvg, 'filled')
axis image
figure;scatter3(coor(:,1),coor(:,2),coor(:,3), 100, nodeAvgStd, 'filled')
axis image

for i = 1:37
    for j = 1:100
        R = rich_club_bu(scrambled_nets(:,:,i,j) > 0);
        [~,p] = findpeaks(R);
if ~isempty(p)
    optimalRC(i, j) = p(length(p));
else
    optimalRC(i, j) = nan;
end
    end
end
choice = nanmean(nanmean(optimalRC));
STD = mean(nanstd(optimalRC, [], 2));
degree = squeeze(mean(sum(scrambled_nets(:,:,:,1) > 0),3));
keep = degree > choice;
figure; boxplot(nodeAvgAvg, keep)
figure; boxplot(nodeAvgStd, keep)
x = nodeAvgStd(keep);
y = nodeAvgStd(~keep);
results(1) = ttest2(x,y); 
x = nodeAvgAvg(keep);
y = nodeAvgAvg(~keep);
results(2) = ttest2(x,y); 


keep = degree > choice + STD;
figure; boxplot(nodeAvgAvg, keep)
figure; boxplot(nodeAvgStd, keep)
x = nodeAvgStd(keep);
y = nodeAvgStd(~keep);
results(3) = ttest2(x,y); 
x = nodeAvgAvg(keep);
y = nodeAvgAvg(~keep);
results(4) = ttest2(x,y); 

keep = degree > choice - STD;
figure; boxplot(nodeAvgAvg, keep)
figure; boxplot(nodeAvgStd, keep)
x = nodeAvgStd(keep);
y = nodeAvgStd(~keep);
results(5) = ttest2(x,y); 
x = nodeAvgAvg(keep);
y = nodeAvgAvg(~keep);
results(6) = ttest2(x,y); 