%% data registeration
d1 = load('E:\daneshgah\arshad\terme 2 - spring 02\CognitiveScience\final project\data\sub-01_Task_paradigm\sub-01_data_B1');
d2 = load('E:\daneshgah\arshad\terme 2 - spring 02\CognitiveScience\final project\data\sub-01_Task_paradigm\sub-01_data_B2');
d3 = load('E:\daneshgah\arshad\terme 2 - spring 02\CognitiveScience\final project\data\sub-01_Task_paradigm\sub-01_data_B3');
d4 = load('E:\daneshgah\arshad\terme 2 - spring 02\CognitiveScience\final project\data\sub-01_Task_paradigm\sub-01_data_B4');
%% extracting match nonmatch indices for each category of stimuli
%%% B1conds
d11 =d1.tmpB.trialSet;
d111 = extractfield(d11, 'simgName');
d1trial = extractfield(d11, 'trialType');
face = ["f", "I"];
match = "Match"; 
B1facematch = transpose(find(startsWith(d111, face)== 1 & strcmp(d1trial, match) == 1));
B1facenonmatch = transpose(find(startsWith(d111, face)== 1 & strcmp(d1trial, match) == 0));
B1nonfacematch = transpose(find(startsWith(d111, face)== 0 & strcmp(d1trial, match) == 1));
B1nonfacenonmatch = transpose(find(startsWith(d111, face)== 0 & strcmp(d1trial, match) == 0));
%% B2conds
d22 =d2.tmpB.trialSet;
d222 = extractfield(d22, 'simgName');
d2trial = extractfield(d22, 'trialType');
face = ["f", "I"];
match = "Match"; 
B2facematch = transpose(find(startsWith(d222, face)== 1 & strcmp(d2trial, match) == 1));
B2facenonmatch = transpose(find(startsWith(d222, face)== 1 & strcmp(d2trial, match) == 0));
B2nonfacematch = transpose(find(startsWith(d222, face)== 0 & strcmp(d2trial, match) == 1));
B2nonfacenonmatch = transpose(find(startsWith(d222, face)== 0 & strcmp(d2trial, match) == 0));
%% B3conditions
d33 =d3.tmpB.trialSet;
d333 = extractfield(d33, 'simgName');
face = ["f", "I"];
d3trial = extractfield(d33, 'trialType');
match = "Match"; 
B3facematch = transpose(find(startsWith(d333, face)== 1 & strcmp(d3trial, match) == 1));
B3facenonmatch = transpose(find(startsWith(d333, face)== 1 & strcmp(d3trial, match) == 0));
B3nonfacematch = transpose(find(startsWith(d333, face)== 0 & strcmp(d3trial, match) == 1));
B3nonfacenonmatch = transpose(find(startsWith(d333, face)== 0 & strcmp(d3trial, match) == 0));
%% B4conditions
d44 =d4.tmpB.trialSet;
d444 = extractfield(d44, 'simgName');
face = ["f", "I"];
d4trial = extractfield(d44, 'trialType');
match = "Match"; 
B4facematch = transpose(find(startsWith(d444, face)== 1 & strcmp(d4trial, match) == 1));
B4facenonmatch = transpose(find(startsWith(d444, face)== 1 & strcmp(d4trial, match) == 0));
B4nonfacematch = transpose(find(startsWith(d444, face)== 0 & strcmp(d4trial, match) == 1));
B4nonfacenonmatch = transpose(find(startsWith(d444, face)== 0 & strcmp(d4trial, match) == 0));
%% FFA ttest 1
[h1, p1, ci1, stats1] = ttest(BetaNonmatchf, Betamatchf);
nu1 = stats1.df; tval1 = stats1.tstat; tcrit1 = tinv(0.95, nu1); tvalpdf1 = tpdf(tval1, nu1);
k = linspace(-20,20,600);
tdistpdf1 = tpdf(k, nu1);
%% 2
mean1 = mean(BetaNonmatchf) -  mean(Betamatchf);
histfit(BetaNonmatchf - Betamatchf);
title('histogram of Beta weights of the contrast matching- nonmatching faces with a fitted t-dist ')
%% plotting/ 3
figure
plot(k, tdistpdf1)
hold on
scatter(tval1, tvalpdf1, "filled")
x_points1 = [mean1-ci1(1), mean1-ci1(1), ci1(2)+mean1, ci1(2)+mean1];
y_points1 = [0, 0.4, 0.4, 0];
xline(tcrit1, "--");
xline(mean1, 'r');
 fill(x_points1, y_points1, [0 0.4470 0.7410], 'FaceAlpha', 0.1, 'EdgeColor', 'none') 
legend(["NonMatch - Match t pdf", "t-statistic", "Critical t Cutoff",...
    "difference in mean for the two vectors", "CI of mean differences"])
xlim([-5, 7])
title("t-dist of the difference in BetaWeights with the tstats illustrated")
%% PPA ttest /1
[h2, p2, ci2, stats2] = ttest(BetaNonMatchobj, BetaMatchobj);
nu2 = stats2.df; tval2 = stats2.tstat; tcrit2 = tinv(0.95, nu2); tvalpdf2 = tpdf(tval2, nu2);
%k = linspace(-20,20,600);
tdistpdf2 = tpdf(k, nu2);
%% 2
mean2 = mean(BetaNonMatchobj) -  mean(BetaMatchobj);
histfit(BetaNonMatchobj-  BetaMatchobj);
title('histogram of Beta weights of the contrast matching- nonmatching NonFaces with a fitted t-dist ')
%% plotting/ 3
figure
plot(k, tdistpdf2)
hold on
scatter(tval2, tvalpdf2, "filled")
x_points2 = [mean2-ci2(1), mean2-ci2(1), ci2(2)+mean2, ci2(2)+mean2];
y_points2 = [0, 0.4, 0.4, 0];
xline(tcrit2, "--");
xline(mean2, 'r');
 fill(x_points2, y_points2, [0 0.4470 0.7410], 'FaceAlpha', 0.1, 'EdgeColor', 'none') 
legend(["NonMatch - Match t pdf", "t-statistic", "Critical t Cutoff",...
    "difference in mean for the two vectors", "CI of mean differences"])
xlim([-5, 7])
title("t-dist of the difference in BetaWeights with the tstats illustrated")
%% PPA vs FFA
face = BetaNonmatchf - Betamatchf;
nonface = BetaNonMatchobj - BetaMatchobj;
[h3, p3, ci3, stats3] = ttest(face,nonface, 'Tail', 'Right');
nu3 = stats3.df; tval3 = stats3.tstat; tcrit3 = tinv(0.95, nu3); tvalpdf3 = tpdf(tval3, nu3);
%k = linspace(-20,20,600);
tdistpdf3 = tpdf(k, nu3);
%% 2
mean3 = mean(face) -  mean(nonface);
histfit(face-nonface);
title('histogram of Beta weights vectors obtained from "nonmatch-match" trials for the face - nonface categories with a fitted t-dist ')
%% plotting/ 3
figure
plot(k, tdistpdf3)
hold on
scatter(tval3, tvalpdf3, "filled")
x_points3 = [mean3-ci3(1), mean3-ci3(1), 10,  10];
y_points3 = [0, 0.4, 0.4, 0];
xline(tcrit3, "--");
xline(mean3, 'r');
 fill(x_points3, y_points3, [0 0.4470 0.7410], 'FaceAlpha', 0.1, 'EdgeColor', 'none') 
legend(["face - nonface t pdf", "t-statistic", "Critical t Cutoff",...
    "difference in mean for the two vectors", "CI of mean differences"])
xlim([-5, 7])
title("t-dist of the difference in BetaWeights with the tstats illustrated")
