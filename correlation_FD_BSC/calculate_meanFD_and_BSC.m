
%AAL Atlas designations (from
%http://doc.pmod.com/pneuro/pneuro.html?aal-mergedatlas6752.html) 

frontal=[1:28,69,70];
insula_cing=[29:36];
temporal=[39,40,79:90];
occipital=[43:54];
parietal=[55:68];

L_frontal=frontal(logical(mod(frontal,2)))
R_frontal=frontal(logical(~mod(frontal,2)))

L_insula_cing=insula_cing(logical(mod(insula_cing,2)))
R_insula_cing=insula_cing(logical(~mod(insula_cing,2)))

L_temporal=temporal(logical(mod(temporal,2)))
R_temporal=temporal(logical(~mod(temporal,2)))

L_occipital=occipital(logical(mod(occipital,2)))
R_occipital=occipital(logical(~mod(occipital,2)))

L_parietal=parietal(logical(mod(parietal,2)))
R_parietal=parietal(logical(~mod(parietal,2)))

lnames=["L_frontal","L_insula_cing","L_temporal","L_occipital","L_parietal"];
rnames=["R_frontal","R_insula_cing","R_temporal","R_occipital","R_parietal"]
list=[frontal,insula_cing,temporal,occipital,parietal];
llist=[{L_frontal},{L_insula_cing},{L_temporal},{L_occipital},{L_parietal}];
rlist=[{R_frontal},{R_insula_cing},{R_temporal},{R_occipital},{R_parietal}];

%load AAL labels
AAL_L=load('/Users/emilyolafson/Documents/BSC/AAL_atlas_left.txt'); %40962 with atlas labels for each vertex.
AAL_L=AAL_L(1:40962);
AAL_R=load('/Users/emilyolafson/Documents/BSC/AAL_atlas_right.txt'); %40962 with atlas labels for each vertex.
AAL_R=AAL_R(1:40962);
% from: /opt/quarantine/CIVET/1.1.12/build/CIVET-1.1.12/models/

%% compute mean BSC and FD for each site.

%KKI subject ids
kki=csvread('BSC/kki.csv')

%FD
for i=1:size(kki,1)
    params = load(strcat('/Users/emilyolafson/Documents/BSC/motion_params/KKI/QC_FDpower_rp_rest.nii.gz_', num2str(kki(i)), '.mat'));
    params=params.R;
    a=mean(params);
    fid=fopen(strcat('/Users/emilyolafson/Documents/BSC/motion_params/KKI/',num2str(kki(i)), '_meanFD.txt'), 'w');
    fprintf(fid,'%f', a)
end

%BSC (global & lobar)
for i=1:size(kki,1)
    bsc=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/KKI/ABIDE_KKI_', num2str(kki(i)), '_model_c_resid_right.txt'));
    a=mean(bsc);
    fid=fopen(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/KKI/ABIDE_KKI_',num2str(kki(i)), '_meanbsc_right.txt'), 'w');
    fprintf(fid,'%f', a)
    for j=1:5
        lobe_vertices=AAL_R==rlist{j};
        lobe_vertices=sum(lobe_vertices,2);
        mean_bsc_lobar_right{i,j}=mean(bsc(logical(lobe_vertices)));
                a=mean(bsc(logical(lobe_vertices)));

        fid=fopen(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/KKI/ABIDE_KKI_',num2str(kki(i)), '_meanbsc', rnames(j), '_right.txt'), 'w');
        fprintf(fid,'%f', a)
        clear lobe_vertices
    end
end

for i=1:size(kki,1)
    bsc=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/KKI/ABIDE_KKI_', num2str(kki(i)), '_model_c_resid_left.txt'));
    a=mean(bsc);
    fid=fopen(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/KKI/ABIDE_KKI_',num2str(kki(i)), '_meanbsc_left.txt'), 'w');
    fprintf(fid,'%f', a)
    for j=1:5
        lobe_vertices=AAL_L==llist{j};
        lobe_vertices=sum(lobe_vertices,2);
        mean_bsc_lobar_left{i,j}=mean(bsc(logical(lobe_vertices)));
                a=mean(bsc(logical(lobe_vertices)));

        fid=fopen(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/KKI/ABIDE_KKI_',num2str(kki(i)), '_meanbsc', rnames(j), '_left.txt'), 'w');
        fprintf(fid,'%f', a)
        clear lobe_vertices
    end
end

mean_bsc_lobar_right;
mean_bsc_lobar_left;

% OHSU subject ids
%ohsu=csvread('BSC/ohsu.csv')

for i=1:size(ohsu,1)
    params = load(strcat('/Users/emilyolafson/Documents/BSC/motion_params/OHSU/QC_FDpower_rp_1_', num2str(ohsu(i)), '_rest.mat'));
    params=params.R;
    a=mean(params);
    fid=fopen(strcat('/Users/emilyolafson/Documents/BSC/motion_params/OHSU/',num2str(ohsu(i)), '_meanFD.txt'), 'w');
    fprintf(fid,'%f', a)
end

%BSC (global & lobar)
for i=1:size(ohsu,1)
    bsc=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/OHSU/ABIDEII_OHSU_1_', num2str(ohsu(i)), '_model_c_resid_right.txt'));
    a=mean(bsc);
    fid=fopen(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/OHSU/ABIDEII_OHSU_1_',num2str(ohsu(i)), '_meanbsc_right.txt'), 'w');
    fprintf(fid,'%f', a)
    fclose(fid)
    for j=1:5
        lobe_vertices=AAL_R==rlist{j};
        lobe_vertices=sum(lobe_vertices,2);
        mean_bsc_lobar_right{i,j}=mean(bsc(logical(lobe_vertices)));
        a=mean(bsc(logical(lobe_vertices)));
        fid=fopen(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/OHSU/ABIDEII_OHSU_1_',num2str(ohsu(i)), '_meanbsc', rnames(j), '_right.txt'), 'w');
        fprintf(fid,'%f', a)
        fclose(fid)
        clear lobe_vertices
    end
end

for i=1:size(ohsu,1)
    bsc=load(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/OHSU/ABIDEII_OHSU_1_', num2str(ohsu(i)), '_model_c_resid_left.txt'));
    a=mean(bsc);
    fid=fopen(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/OHSU/ABIDEII_OHSU_1_',num2str(ohsu(i)), '_meanbsc_left.txt'), 'w');
    fprintf(fid,'%f', a)
        fclose(fid)
    for j=1:5
        lobe_vertices=AAL_L==llist{j};
        lobe_vertices=sum(lobe_vertices,2);
        mean_bsc_lobar_left{i,j}=mean(bsc(logical(lobe_vertices)));
        a=mean(bsc(logical(lobe_vertices)));
        fid=fopen(strcat('/Users/emilyolafson/Documents/BSC/bsc_values/OHSU/ABIDEII_OHSU_1_',num2str(ohsu(i)), '_meanbsc', rnames(j), '_left.txt'), 'w');
        fprintf(fid,'%f', a)
        fclose(fid)
        clear lobe_vertices
    end
end