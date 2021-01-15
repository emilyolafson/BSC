% weighted AIC 

%load site size weights
gf=readtable('gf.csv')

kki=sum(contains(gf.Site_combined, 'KKI'))
nyu=sum(contains(gf.Site_combined, 'NYU'))
ohsu=sum(contains(gf.Site_combined, 'OHSU'))
to=sum(contains(gf.Site_combined, 'TORONTO'))
sdsu=sum(contains(gf.Site_combined, 'SDSU'))
cam=sum(contains(gf.Site_combined, 'CAM'))
cambridge=sum(contains(gf.Site_combined, 'Cambridge'))
um=sum(contains(gf.Site_combined, 'UM'))
maxmun=sum(contains(gf.Site_combined, 'MAX_MUN'))
IP=sum(contains(gf.Site_combined, 'IP'))
iop=sum(contains(gf.Site_combined, 'IoP'))

weights=[kki,nyu,ohsu,to,sdsu,cam,cambridge,um,maxmun,IP,iop];
weights=weights./sum(weights)

% load aic results
aic_sex_L=readtable('AIC_sex/sex_whichbest_L_combined_model_c.csv');
aic_sex_L=table2array(aic_sex_L);
aic_sex_R=readtable('AIC_sex/sex_whichbest_R_combined_model_c.csv');
aic_sex_R=table2array(aic_sex_R);
aic_age_L=readtable('AIC_age/age_whichbest_L_combined_model_c.csv');
aic_age_L=table2array(aic_age_L);
aic_age_R=readtable('AIC_age/age_whichbest_R_combined_model_c.csv');
aic_age_R=table2array(aic_age_R);

%SEX
model1=aic_sex_L==1;
model1=model1*weights'
writematrix(model1, '/Users/emilyolafson/Documents/BSC/aic_noitxn/sex_L_weighted_model1.txt')
model1=aic_sex_R==1;
model1=model1*weights'
writematrix(model1, '/Users/emilyolafson/Documents/BSC/aic_noitxn/sex_R_weighted_model1.txt')

model2=aic_sex_L==2;
model2=model2*weights';
writematrix(model2, '/Users/emilyolafson/Documents/BSC/aic_noitxn/sex_L_weighted_model2.txt')
model2=aic_sex_R==2;
model2=model2*weights';
writematrix(model2, '/Users/emilyolafson/Documents/BSC/aic_noitxn/sex_R_weighted_model2.txt')

model3=aic_sex_L==3;
model3=model3*weights'
writematrix(model3, '/Users/emilyolafson/Documents/BSC/aic_noitxn/sex_L_weighted_model3.txt')
model3=aic_sex_R==3;
model3=model3*weights'
writematrix(model3, '/Users/emilyolafson/Documents/BSC/aic_noitxn/sex_R_weighted_model3.txt')



%AGE
model1=aic_age_L==1;
model1=model1*weights'
writematrix(model1, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age_L_weighted_model1.txt')
model1=aic_age_R==1;
model1=model1*weights'
writematrix(model1, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age_R_weighted_model1.txt')

model2=aic_age_L==2;
model2=model2*weights';
writematrix(model2, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age_L_weighted_model2.txt')
model2=aic_age_R==2;
model2=model2*weights';
writematrix(model2, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age_R_weighted_model2.txt')

model3=aic_age_L==3;
model3=model3*weights'
writematrix(model3, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age_L_weighted_model3.txt')
model3=aic_age_R==3;
model3=model3*weights'
writematrix(model3, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age_R_weighted_model3.txt')



kki=aic_sex_L(:,1);
writematrix(kki, '/Users/emilyolafson/Documents/BSC/aic_noitxn/kki_L.txt')

nyu=aic_sex_L(:,2);
writematrix(nyu, '/Users/emilyolafson/Documents/BSC/aic_noitxn/nyu_L.txt')

ohsu=aic_sex_L(:,3);
writematrix(ohsu, '/Users/emilyolafson/Documents/BSC/aic_noitxn/ohsu_L.txt')

to=aic_sex_L(:,4);
writematrix(to, '/Users/emilyolafson/Documents/BSC/aic_noitxn/to_L.txt')

sdsu=aic_sex_L(:,5);
writematrix(sdsu, '/Users/emilyolafson/Documents/BSC/aic_noitxn/sdsu_L.txt')

cam=aic_sex_L(:,6);
writematrix(cam, '/Users/emilyolafson/Documents/BSC/aic_noitxn/cam_L.txt')

cambridge=aic_sex_L(:,7);
writematrix(cambridge, '/Users/emilyolafson/Documents/BSC/aic_noitxn/cambridge_L.txt')

um=aic_sex_L(:,8);
writematrix(um, '/Users/emilyolafson/Documents/BSC/aic_noitxn/um_L.txt')

maxmun=aic_sex_L(:,9);
writematrix(maxmun, '/Users/emilyolafson/Documents/BSC/aic_noitxn/maxmun_L.txt')

IP=aic_sex_L(:,10);
writematrix(IP, '/Users/emilyolafson/Documents/BSC/aic_noitxn/IP_L.txt')

iop=aic_sex_L(:,11);
writematrix(iop, '/Users/emilyolafson/Documents/BSC/aic_noitxn/iop_L.txt')


kki=aic_sex_R(:,1);
writematrix(kki, '/Users/emilyolafson/Documents/BSC/aic_noitxn/kki_R.txt')

nyu=aic_sex_R(:,2);
writematrix(nyu, '/Users/emilyolafson/Documents/BSC/aic_noitxn/nyu_R.txt')

ohsu=aic_sex_R(:,3);
writematrix(ohsu, '/Users/emilyolafson/Documents/BSC/aic_noitxn/ohsu_R.txt')

to=aic_sex_R(:,4);
writematrix(to, '/Users/emilyolafson/Documents/BSC/aic_noitxn/to_R.txt')

sdsu=aic_sex_R(:,5);
writematrix(sdsu, '/Users/emilyolafson/Documents/BSC/aic_noitxn/sdsu_R.txt')

cam=aic_sex_R(:,6);
writematrix(cam, '/Users/emilyolafson/Documents/BSC/aic_noitxn/cam_R.txt')

cambridge=aic_sex_R(:,7);
writematrix(cambridge, '/Users/emilyolafson/Documents/BSC/aic_noitxn/cambridge_R.txt')

um=aic_sex_R(:,8);
writematrix(um, '/Users/emilyolafson/Documents/BSC/aic_noitxn/um_R.txt')

maxmun=aic_sex_R(:,9);
writematrix(maxmun, '/Users/emilyolafson/Documents/BSC/aic_noitxn/maxmun_R.txt')

IP=aic_sex_R(:,10);
writematrix(IP, '/Users/emilyolafson/Documents/BSC/aic_noitxn/IP_R.txt')

iop=aic_sex_R(:,11);
writematrix(iop, '/Users/emilyolafson/Documents/BSC/aic_noitxn/iop_R.txt')





kki=aic_age_L(:,1);
writematrix(kki, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/kki_L.txt')

nyu=aic_age_L(:,2);
writematrix(nyu, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/nyu_L.txt')

ohsu=aic_age_L(:,3);
writematrix(ohsu, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/ohsu_L.txt')

to=aic_age_L(:,4);
writematrix(to, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/to_L.txt')

sdsu=aic_age_L(:,5);
writematrix(sdsu, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/sdsu_L.txt')

cam=aic_age_L(:,6);
writematrix(cam, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/cam_L.txt')

cambridge=aic_age_L(:,7);
writematrix(cambridge, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/cambridge_L.txt')

um=aic_age_L(:,8);
writematrix(um, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/um_L.txt')

maxmun=aic_age_L(:,9);
writematrix(maxmun, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/maxmun_L.txt')

IP=aic_age_L(:,10);
writematrix(IP, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/IP_L.txt')

iop=aic_age_L(:,11);
writematrix(iop, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/iop_L.txt')


kki=aic_age_R(:,1);
writematrix(kki, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/kki_R.txt')

nyu=aic_age_R(:,2);
writematrix(nyu, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/nyu_R.txt')

ohsu=aic_age_R(:,3);
writematrix(ohsu, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/ohsu_R.txt')

to=aic_age_R(:,4);
writematrix(to, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/to_R.txt')

sdsu=aic_age_R(:,5);
writematrix(sdsu, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/sdsu_R.txt')

cam=aic_age_R(:,6);
writematrix(cam, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/cam_R.txt')

cambridge=aic_age_R(:,7);
writematrix(cambridge, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/cambridge_R.txt')

um=aic_age_R(:,8);
writematrix(um, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/um_R.txt')

maxmun=aic_age_R(:,9);
writematrix(maxmun, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/maxmun_R.txt')

IP=aic_age_R(:,10);
writematrix(IP, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/IP_R.txt')

iop=aic_age_R(:,11);
writematrix(iop, '/Users/emilyolafson/Documents/BSC/aic_noitxn/age/iop_R.txt')
