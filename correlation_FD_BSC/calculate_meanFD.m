kki=csvread('BSC/kki.csv')

for i=1:size(kki,1)
    params = load(strcat('/Users/emilyolafson/Documents/BSC/motion_params/KKI/QC_FDpower_rp_rest.nii.gz_', num2str(kki(i)), '.mat'));
    params=params.R;
    a=mean(params);
    fid=fopen(strcat('/Users/emilyolafson/Documents/BSC/motion_params/KKI/',num2str(kki(i)), '_meanFD.txt'), 'w');
    fprintf(fid,'%f', a)
end

ohsu=csvread('BSC/ohsu.csv')

for i=1:size(ohsu,1)
    params = load(strcat('/Users/emilyolafson/Documents/BSC/motion_params/OHSU/QC_FDpower_rp_1_', num2str(ohsu(i)), '_rest.mat'));
    params=params.R;
    a=mean(params);
    fid=fopen(strcat('/Users/emilyolafson/Documents/BSC/motion_params/OHSU/',num2str(ohsu(i)), '_meanFD.txt'), 'w');
    fprintf(fid,'%f', a)
end

