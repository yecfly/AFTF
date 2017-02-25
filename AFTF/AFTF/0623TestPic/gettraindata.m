cd('D:\YECWorkspace\CODES\AFTF\AFTF\0623TestPic');
% G1=load('D:\YECWorkspace\CODES\AFTF\AFTF\0623TestPic\G1_20160624_1139_Ransac4\_summary');
% G2=load('D:\YECWorkspace\CODES\AFTF\AFTF\0623TestPic\G2_20160624_1140_Ransac4\_summary');
% G3=load('D:\YECWorkspace\CODES\AFTF\AFTF\0623TestPic\G3_20160624_1141_Ransac4\_summary');
% G4=load('D:\YECWorkspace\CODES\AFTF\AFTF\0623TestPic\G4_20160624_1143_Ransac4\_summary');
% save('G1-G4.mat','G1','G2','G3','G4');
load('G1-G4.mat');

rotation1=G1(:,3)-G1(:,2);
x1=G1(:,5)-G1(:,4);
y1=G1(:,7)-G1(:,6);
ree1=G1(:,9);
sinv1=G1(:,20);
cosv1=G1(:,19);
scale1=G1(:,21);
cf1=G1(:,13);
vx1=G1(:,10);
vy1=G1(:,11);
convin1=G1(:,12);
ax1=G1(:,17);
ay1=G1(:,18);
c1=[G1(:,1)';G1(:,22)';G1(:,23)']';

rotation2=G2(:,3)-G2(:,2);
x2=G2(:,5)-G2(:,4);
y2=G2(:,7)-G2(:,6);
ree2=G2(:,9);
sinv2=G2(:,20);
cosv2=G2(:,19);
scale2=G2(:,21);
cf2=G2(:,13);
vx2=G2(:,10);
vy2=G2(:,11);
convin2=G2(:,12);
ax2=G2(:,17);
ay2=G2(:,18);
c2=[G2(:,1)';G2(:,22)';G2(:,23)']';

rotation3=G3(:,3)-G3(:,2);
x3=G3(:,5)-G3(:,4);
y3=G3(:,7)-G3(:,6);
ree3=G3(:,9);
sinv3=G3(:,20);
cosv3=G3(:,19);
scale3=G3(:,21);
cf3=G3(:,13);
vx3=G3(:,10);
vy3=G3(:,11);
convin3=G3(:,12);
ax3=G3(:,17);
ay3=G3(:,18);
c3=[G3(:,1)';G3(:,22)';G3(:,23)']';

rotation4=G4(:,3)-G4(:,2);
x4=G4(:,5)-G4(:,4);
y4=G4(:,7)-G4(:,6);
ree4=G4(:,9);
sinv4=G4(:,20);
cosv4=G4(:,19);
scale4=G4(:,21);
cf4=G4(:,13);
vx4=G4(:,10);
vy4=G4(:,11);
convin4=G4(:,12);
ax4=G4(:,17);
ay4=G4(:,18);
c4=[G4(:,1)';G4(:,22)';G4(:,23)']';

rotation=[rotation1;rotation2;rotation3;rotation4];
x=[x1;x2;x3;x4];
y=[y1;y2;y3;y4];
ree=[ree1;ree2;ree3;ree4];
sinv=[sinv1;sinv2;sinv3;sinv4];
cosv=[cosv1;cosv2;cosv3;cosv4];
scale=[scale1;scale2;scale3;scale4];
cf=[cf1;cf2;cf3;cf4];
vx=[vx1;vx2;vx3;vx4];
vy=[vy1;vy2;vy3;vy4];
convin=[convin1;convin2;convin3;convin4];
ax=[ax1;ax2;ax3;ax4];
ay=[ay1;ay2;ay3;ay4];
c=[c1;c2;c3;c4];
r_vs_cf=[ree';cf']';

%%
clear rotation1 rotation2 rotation3 rotation4 x1 x2 x3 x4 y1 y2 y3 y4
clear ree1 ree2 ree3 ree4 sinv1 sinv2 sinv3 sinv4 cosv1 cosv2 cosv3 cosv4 
clear scale1 scale2 scale3 scale4 cf1 cf2 cf3 cf4 vx1 vx2 vx3 vx4
clear vy1 vy2 vy3 vy4 convin1 convin2 convin3 convin4 ax1 ax2 ax3 ax4
clear ay1 ay2 ay3 ay4 c1 c2 c3 c4