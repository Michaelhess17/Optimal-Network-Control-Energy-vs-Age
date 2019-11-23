load('C:\Users\Micha\Desktop\Computer Science\MATLAB\nki+controllability\output\nodeAvgResAD.mat')
load('C:\Users\Micha\Desktop\Computer Science\MATLAB\nki+controllability\mat\nki_mats.mat', 'coor', 'ci', 'SC')
nodeAvgResAD = reshape(nodeAvgResAD,[214,663,81]);
nodeAvgAvg = nanmean(nodeAvgResAD, 3);
nodeAvgAvg = nanmean(nodeAvgAvg, 2);
nodeAvgStd = nanmean(std(nodeAvgResAD, [], 3),2);
figure;scatter3(coor(:,1),coor(:,2),coor(:,3), 100, nodeAvgAvg, 'filled')
axis image
figure;scatter3(coor(:,1),coor(:,2),coor(:,3), 100, nodeAvgStd, 'filled')
axis image

for i = 1:663
R = rich_club_bu(SC(:,:,i) > 0);
 [~,p] = findpeaks(R);
if ~isempty(p)
    optimalRC(i) = p(length(p));
else
    optimalRC(i) = nan;
end
end
choice = nanmean(optimalRC);
degree = squeeze(mean(sum(SC > 0,2),3));
STD = nanstd(optimalRC, [], 2);
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