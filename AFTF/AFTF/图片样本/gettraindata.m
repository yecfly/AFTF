% HR=load('I:\DATASETS\¹ãÖÝ°¬²© Í¼Æ¬Ñù±¾\µÚÒ»×é\T5HomoRIGI_20160422_1454 summary');

% HR=load('K:\°¬²© Í¼Æ¬²âÊÔ\G1_20160425_2022_Ransac4\_summary Ptrain');
% xd=HR(:,13);
% yd=HR(:,14);
% rotation=HR(:,2)-HR(:,1);
% ax=HR(:,15);
% ay=HR(:,16);
% sinv=HR(:,18);
% cosv=HR(:,17);
% scale=HR(:,19);
% 
% HR=load('K:\°¬²© Í¼Æ¬²âÊÔ\G2_20160425_2024_Ransac4\_summary Ptrain');
% 
% xd1=HR(:,13);
% yd1=HR(:,14);
% rotation1=HR(:,2)-HR(:,1);
% ax1=HR(:,15);
% ay1=HR(:,16);
% sinv1=HR(:,18);
% cosv1=HR(:,17);
% scale1=HR(:,19);


xdt=[xd;xd1];
ydt=[yd;yd1];
rotationt=[rotation;rotation1];
axt=[ax;ax1];
ayt=[ay;ay1];
sinvt=[sinv;sinv1];
cosvt=[cosv;cosv1];
scalet=[scale;scale1];