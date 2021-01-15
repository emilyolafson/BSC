% Relationship between BSC and FD.
clear ohsu_*
ohsu_demog=readtable('BSC/ohsudemog.txt');
ohsu_demog=ohsu_demog(1:74,:);
ohsu_SUBId=ohsu_demog.Var8;
ohsu_age=ohsu_demog.Var12;
ohsu_sex=ohsu_demog.Var13;
ohsu_dx=ohsu_demog.Var10;


kki_demog=readtable('BSC/kkidemog.txt');
kki_SUBId=kki_demog.Var8;
kki_age=kki_demog.Var12;
kki_sex=kki_demog.Var13;
kki_dx=kki_demog.Var10;

kki=csvread('BSC/kki.csv');
ohsu=csvread('BSC/ohsu.csv');
ohsu=ohsu(1:74);

kki=kki_demog.Var2;
ohsu=ohsu_demog.Var8;



cc=load('BSC/CIVET-CC-mask.txt'); % 1's where there is no midline. 

clear kki_FD
clear kki_BSC_left
clear kki_BSC_*

for i=1:size(kki,1)
    kki_FD{i}=load(strcat('/Users/emilyolafson/Documents/BSC/motion_params/KKI/',num2str(kki(i)), '_meanFD.txt'));
    kki_BSC_left{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/KKI/ABIDE_KKI_', num2str(kki(i)), '_meanbsc_left.txt'));
    kki_BSC_right{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/KKI/ABIDE_KKI_', num2str(kki(i)), '_meanbsc_right.txt'));
    kki_BSC_frontal_left{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/KKI/ABIDE_KKI_', num2str(kki(i)), '_meanbscR_frontal_left.txt'));
    kki_BSC_insula_cing_left{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/KKI/ABIDE_KKI_', num2str(kki(i)), '_meanbscR_insula_cing_left.txt'));
    kki_BSC_occipital_left{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/KKI/ABIDE_KKI_', num2str(kki(i)), '_meanbscR_occipital_left.txt'));
    kki_BSC_parietal_left{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/KKI/ABIDE_KKI_', num2str(kki(i)), '_meanbscR_parietal_left.txt'));
    kki_BSC_temporal_left{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/KKI/ABIDE_KKI_', num2str(kki(i)), '_meanbscR_temporal_left.txt'));
    kki_BSC_frontal_right{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/KKI/ABIDE_KKI_', num2str(kki(i)), '_meanbscR_frontal_right.txt'));
    kki_BSC_insula_cing_right{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/KKI/ABIDE_KKI_', num2str(kki(i)), '_meanbscR_insula_cing_right.txt'));
    kki_BSC_occipital_right{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/KKI/ABIDE_KKI_', num2str(kki(i)), '_meanbscR_occipital_right.txt'));
    kki_BSC_parietal_right{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/KKI/ABIDE_KKI_', num2str(kki(i)), '_meanbscR_parietal_right.txt'));
    kki_BSC_temporal_right{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/KKI/ABIDE_KKI_', num2str(kki(i)), '_meanbscR_temporal_right.txt'));
end

clear ohsu_FD
clear ohsu_BSC*

for i=1:200
    ohsu_FD{i}=load(strcat('/Users/emilyolafson/Documents/BSC/motion_params/OHSU/',num2str(ohsu(i)), '_meanFD.txt'));
    ohsu_BSC_left{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/OHSU/ABIDEII_OHSU_1_', num2str(ohsu(i)), '_meanbsc_left.txt'));
    ohsu_BSC_right{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/OHSU/ABIDEII_OHSU_1_', num2str(ohsu(i)), '_meanbsc_right.txt'));
    ohsu_BSC_frontal_left{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/OHSU/ABIDEII_OHSU_1_', num2str(ohsu(i)), '_meanbscR_frontal_left.txt'));
    ohsu_BSC_insula_cing_left{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/OHSU/ABIDEII_OHSU_1_', num2str(ohsu(i)), '_meanbscR_insula_cing_left.txt'));
    ohsu_BSC_occipital_left{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/OHSU/ABIDEII_OHSU_1_', num2str(ohsu(i)), '_meanbscR_occipital_left.txt'));
    ohsu_BSC_parietal_left{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/OHSU/ABIDEII_OHSU_1_', num2str(ohsu(i)), '_meanbscR_parietal_left.txt'));
    ohsu_BSC_temporal_left{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/OHSU/ABIDEII_OHSU_1_', num2str(ohsu(i)), '_meanbscR_temporal_left.txt'));
    ohsu_BSC_frontal_right{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/OHSU/ABIDEII_OHSU_1_', num2str(ohsu(i)), '_meanbscR_frontal_right.txt'));
    ohsu_BSC_insula_cing_right{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/OHSU/ABIDEII_OHSU_1_', num2str(ohsu(i)), '_meanbscR_insula_cing_right.txt'));
    ohsu_BSC_occipital_right{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/OHSU/ABIDEII_OHSU_1_', num2str(ohsu(i)), '_meanbscR_occipital_right.txt'));
    ohsu_BSC_parietal_right{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/OHSU/ABIDEII_OHSU_1_', num2str(ohsu(i)), '_meanbscR_parietal_right.txt'));
    ohsu_BSC_temporal_right{i}=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/OHSU/ABIDEII_OHSU_1_', num2str(ohsu(i)), '_meanbscR_temporal_right.txt'));
end


kki_FD=cell2mat(kki_FD)
kki_BSC_left=cell2mat(kki_BSC_left)
kki_BSC_right=cell2mat(kki_BSC_right)
kki_BSC_frontal_left=cell2mat(kki_BSC_frontal_left)
kki_BSC_frontal_right=cell2mat(kki_BSC_frontal_right)
kki_BSC_insula_cing_left=cell2mat(kki_BSC_insula_cing_left)
kki_BSC_insula_cing_right=cell2mat(kki_BSC_insula_cing_right)
kki_BSC_occipital_left=cell2mat(kki_BSC_occipital_left)
kki_BSC_occipital_right=cell2mat(kki_BSC_occipital_right)
kki_BSC_parietal_left=cell2mat(kki_BSC_parietal_left)
kki_BSC_parietal_right=cell2mat(kki_BSC_parietal_right)
kki_BSC_temporal_left=cell2mat(kki_BSC_temporal_left)
kki_BSC_temporal_right=cell2mat(kki_BSC_temporal_right)


ohsu_FD=cell2mat(ohsu_FD)
ohsu_BSC_left=cell2mat(ohsu_BSC_left)
ohsu_BSC_right=cell2mat(ohsu_BSC_right)
ohsu_BSC_frontal_left=cell2mat(ohsu_BSC_frontal_left)
ohsu_BSC_frontal_right=cell2mat(ohsu_BSC_frontal_right)
ohsu_BSC_insula_cing_left=cell2mat(ohsu_BSC_insula_cing_left)
ohsu_BSC_insula_cing_right=cell2mat(ohsu_BSC_insula_cing_right)
ohsu_BSC_occipital_left=cell2mat(ohsu_BSC_occipital_left)
ohsu_BSC_occipital_right=cell2mat(ohsu_BSC_occipital_right)
ohsu_BSC_parietal_left=cell2mat(ohsu_BSC_parietal_left)
ohsu_BSC_parietal_right=cell2mat(ohsu_BSC_parietal_right)
ohsu_BSC_temporal_left=cell2mat(ohsu_BSC_temporal_left)
ohsu_BSC_temporal_right=cell2mat(ohsu_BSC_temporal_right)



X=[ohsu_age,ohsu_sex,ohsu_FD'];
y=[ohsu_BSC_frontal_left+ohsu_BSC_frontal_right]./2;
mdl_front=fitlm(X,y);

y=[ohsu_BSC_insula_cing_left+ohsu_BSC_insula_cing_right]./2;
mdl_insula=fitlm(X,y);

y=[ohsu_BSC_occipital_left+ohsu_BSC_occipital_right]./2;
mdl_occip=fitlm(X,y);

y=[ohsu_BSC_parietal_left+ohsu_BSC_parietal_right]./2;
mdl_parietal=fitlm(X,y);

y=[ohsu_BSC_temporal_left+ohsu_BSC_temporal_right]./2;
mdl_temporal=fitlm(X,y);

y=[ohsu_BSC_left+ohsu_BSC_right]./2;
mdl=fitlm(X,y)


X=[kki_age,kki_sex,kki_FD'];
y=[kki_BSC_frontal_left+kki_BSC_frontal_right]./2;
mdl_front=fitlm(X,y);

y=[kki_BSC_insula_cing_left+kki_BSC_insula_cing_right]./2;
mdl_insula=fitlm(X,y);

y=[kki_BSC_occipital_left+kki_BSC_occipital_right]./2;
mdl_occip=fitlm(X,y);

y=[kki_BSC_parietal_left+kki_BSC_parietal_right]./2;
mdl_parietal=fitlm(X,y);

y=[kki_BSC_temporal_left+kki_BSC_temporal_right]./2;
mdl_temporal=fitlm(X,y);

y=[kki_BSC_left+kki_BSC_right]./2;
mdl=fitlm(X,y);

tiledlayout(1,2,'padding', 'none')
nexttile;
plot(ohsu_FD, ohsu_BSC_left, '*b')
title('FD vs. mean left BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')
nexttile;
plot(ohsu_FD, ohsu_BSC_right, '*b')
title('FD vs. mean right BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')
sgtitle('OHSU ABIDE II FD vs. BSC')

tiledlayout(1,2,'padding', 'none')
nexttile;
plot(kki_FD, kki_BSC_left, '*b')
title('FD vs. mean left BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')
nexttile;
plot(kki_FD, kki_BSC_right, '*b')
title('FD vs. mean right BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')
sgtitle('KKI ABIDE I FD vs. BSC')


tiledlayout(5,2, 'padding', 'none')
nexttile;

plot(ohsu_FD, ohsu_BSC_frontal_left, '*b')
title('FD vs. mean left frontal BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')

nexttile;
plot(ohsu_FD, ohsu_BSC_frontal_right, '*b')
title('FD vs. mean right frontal BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')

nexttile;
plot(ohsu_FD, ohsu_BSC_insula_cing_left, '*b')
title('FD vs. mean left insula/cingulate BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')

nexttile;
plot(ohsu_FD, ohsu_BSC_insula_cing_right, '*b')
title('FD vs. mean right insula/cingulate BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')

nexttile;
plot(ohsu_FD, ohsu_BSC_occipital_left, '*b')
title('FD vs. mean left occpital BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')

nexttile;
plot(ohsu_FD, ohsu_BSC_occipital_right, '*b')
title('FD vs. mean right occpital BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')

nexttile;
plot(ohsu_FD, ohsu_BSC_parietal_left, '*b')
title('FD vs. mean left parietal BSC')
xlabel('Average Framewise Displacement')

nexttile;
plot(ohsu_FD, ohsu_BSC_parietal_right, '*b')
title('FD vs. mean right parietal BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')

nexttile;
plot(ohsu_FD, ohsu_BSC_temporal_left, '*b')
title('FD vs. mean left temporal BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')

nexttile;
plot(ohsu_FD, ohsu_BSC_temporal_right, '*b')
title('FD vs. mean right temporal BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')


tiledlayout(5,2, 'padding', 'none')
nexttile;

plot(kki_FD, kki_BSC_frontal_left, '*b')
title('FD vs. mean left frontal BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')

nexttile;
plot(kki_FD, kki_BSC_frontal_right, '*b')
title('FD vs. mean right frontal BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')

nexttile;
plot(kki_FD, kki_BSC_insula_cing_left, '*b')
title('FD vs. mean left insula/cingulate BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')

nexttile;
plot(kki_FD, kki_BSC_insula_cing_right, '*b')
title('FD vs. mean right insula/cingulate BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')

nexttile;
plot(kki_FD, kki_BSC_occipital_left, '*b')
title('FD vs. mean left occpital BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')

nexttile;
plot(kki_FD, kki_BSC_occipital_right, '*b')
title('FD vs. mean right occpital BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')

nexttile;
plot(kki_FD, kki_BSC_parietal_left, '*b')
title('FD vs. mean left parietal BSC')
xlabel('Average Framewise Displacement')

nexttile;
plot(kki_FD, kki_BSC_parietal_right, '*b')
title('FD vs. mean right parietal BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')

nexttile;
plot(kki_FD, kki_BSC_temporal_left, '*b')
title('FD vs. mean left temporal BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')

nexttile;
plot(kki_FD, kki_BSC_temporal_right, '*b')
title('FD vs. mean right temporal BSC')
xlabel('Average Framewise Displacement')
ylabel('Average BSC')
