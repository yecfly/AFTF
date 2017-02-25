% 3
% t=3;
% v1=load('I:\DATASETS\广州艾博 图片样本\Validate一_20160423_1906 summary');
% v2=load('I:\DATASETS\广州艾博 图片样本\Validate二_20160423_1908 summary');
% v3=load('I:\DATASETS\广州艾博 图片样本\Validate三_20160423_1909 summary');
% v4=load('I:\DATASETS\广州艾博 图片样本\Validate四_20160423_1912 summary');

%% 9
% t=9;
% v1=load('I:\DATASETS\广州艾博 图片样本\Validate一_20160423_1921 summary');
% v2=load('I:\DATASETS\广州艾博 图片样本\Validate二_20160423_1922 summary');
% v3=load('I:\DATASETS\广州艾博 图片样本\Validate三_20160423_1923 summary');
% v4=load('I:\DATASETS\广州艾博 图片样本\Validate四_20160423_1926 summary');

%% 1
% t=1;
% v1=load('I:\DATASETS\广州艾博 图片样本\Validate一_20160423_1937 summary');
% v2=load('I:\DATASETS\广州艾博 图片样本\Validate二_20160423_1938 summary');
% v3=load('I:\DATASETS\广州艾博 图片样本\Validate三_20160423_1939 summary');
% v4=load('I:\DATASETS\广州艾博 图片样本\Validate四_20160423_1943 summary');

%%
% v1=load('K:\艾博 图片测试\G1_20160427_2039_Ransac4\_summary');
% v2=load('K:\艾博 图片测试\G2_20160427_2040_Ransac4\_summary');
% v3=load('K:\艾博 图片测试\G3_20160428_1541_Ransac4\_summary');
% v4=load('K:\艾博 图片测试\G4_20160427_2045_Ransac4\_summary');
%% e.rotation poly21 filter7
v1=load('K:\艾博 图片测试\G1_20160428_1538_Ransac4\_summary');
v2=load('K:\艾博 图片测试\G2_20160428_1540_Ransac4\_summary');
v3=load('K:\艾博 图片测试\G3_20160428_1541_Ransac4\_summary');
v4=load('K:\艾博 图片测试\G4_20160428_1544_Ransac4\_summary');
%%
% v1=load('K:\艾博 图片测试\G1_20160428_1559_Ransac4\_summary');
% v2=load('K:\艾博 图片测试\G2_20160428_1601_Ransac4\_summary');
% v3=load('K:\艾博 图片测试\G3_20160428_1602_Ransac4\_summary');
% v4=load('K:\艾博 图片测试\G4_20160428_1605_Ransac4\_summary');
%%
[e1,m,d1]=validSummary(v1);
[e2,m,d2]=validSummary(v2);
[e3,m,d3]=validSummary(v3);
[e4,m,d4]=validSummary(v4);

% rnpm=[d1.rn.pm';d2.rn.pm';d3.rn.pm';d4.rn.pm'];
% rnpv=[d1.rn.pv';d2.rn.pv';d3.rn.pv';d4.rn.pv'];
% rnfp=[d1.rn.fp';d2.rn.fp';d3.rn.fp';d4.rn.fp'];
% rnno=[d1.rn.no';d2.rn.no';d3.rn.no';d4.rn.no'];
% 
% rppm=[d1.rp.pm';d2.rp.pm';d3.rp.pm';d4.rp.pm'];
% rppv=[d1.rp.pv';d2.rp.pv';d3.rp.pv';d4.rp.pv'];
% rpfp=[d1.rp.fp';d2.rp.fp';d3.rp.fp';d4.rp.fp'];
% rpno=[d1.rp.no';d2.rp.no';d3.rp.no';d4.rp.no'];
% 
% xnpm=[d1.xn.pm';d2.xn.pm';d3.xn.pm';d4.xn.pm'];
% xnpv=[d1.xn.pv';d2.xn.pv';d3.xn.pv';d4.xn.pv'];
% xnfp=[d1.xn.fp';d2.xn.fp';d3.xn.fp';d4.xn.fp'];
% xnno=[d1.xn.no';d2.xn.no';d3.xn.no';d4.xn.no'];
% 
% xppm=[d1.xp.pm';d2.xp.pm';d3.xp.pm';d4.xp.pm'];
% xppv=[d1.xp.pv';d2.xp.pv';d3.xp.pv';d4.xp.pv'];
% xpfp=[d1.xp.fp';d2.xp.fp';d3.xp.fp';d4.xp.fp'];
% xpno=[d1.xp.no';d2.xp.no';d3.xp.no';d4.xp.no'];
% 
% % ynpm=[d1.yn.pm';d2.yn.pm';d3.yn.pm';d4.yn.pm'];
% % ynpv=[d1.yn.pv';d2.yn.pv';d3.yn.pv';d4.yn.pv'];
% % ynfp=[d1.yn.fp';d2.yn.fp';d3.yn.fp';d4.yn.fp'];
% % ynno=[d1.yn.no';d2.yn.no';d3.yn.no';d4.yn.no'];
% ynpm=[d3.yn.pm';d4.yn.pm'];
% ynpv=[d3.yn.pv';d4.yn.pv'];
% ynfp=[d3.yn.fp';d4.yn.fp'];
% ynno=[d3.yn.no';d4.yn.no'];
% 
% yppm=[d1.yp.pm';d2.yp.pm';d3.yp.pm';d4.yp.pm'];
% yppv=[d1.yp.pv';d2.yp.pv';d3.yp.pv';d4.yp.pv'];
% ypfp=[d1.yp.fp';d2.yp.fp';d3.yp.fp';d4.yp.fp'];
% ypno=[d1.yp.no';d2.yp.no';d3.yp.no';d4.yp.no'];
% 
% tppm=[d1.tp.pm';d2.tp.pm';d3.tp.pm';d4.tp.pm'];
% tppv=[d1.tp.pv';d2.tp.pv';d3.tp.pv';d4.tp.pv'];
% tpfp=[d1.tp.fp';d2.tp.fp';d3.tp.fp';d4.tp.fp'];
% tpno=[d1.tp.no';d2.tp.no';d3.tp.no';d4.tp.no'];

% figure('Name','Rotation')
% plot3(rnpm,rnpv,rnfp,'bx');
% axis vis3d
% hold on
% plot3(rppm,rppv,rpfp,'g*');
% 
% figure('Name','X')
% plot3(xnpm,xnpv,xnfp,'bx');
% axis vis3d
% hold on
% plot3(xppm,xppv,xpfp,'g*');
% 
% figure('Name','Y')
% plot3(ynpm,ynpv,ynfp,'bx');
% axis vis3d
% hold on
% plot3(yppm,yppv,ypfp,'g*');
% fd=fopen('compare data.txt','a');
% fp'rintf(fd,'Ransac threshold %d 角度误差:%g x误差:%g y误差:%g\n',t,m.rthresh,m.xthresh,m.ythresh);
% fp'rintf(fd,'第1组 %g\t %g\t %g\t %g\t %g\t\n',e1.evapoints,e1.tr,e1.rr,e1.xr,e1.yr);
% fp'rintf(fd,'第2组 %g\t %g\t %g\t %g\t %g\t\n',e2.evapoints,e2.tr,e2.rr,e2.xr,e2.yr);
% fp'rintf(fd,'第3组 %g\t %g\t %g\t %g\t %g\t\n',e3.evapoints,e3.tr,e3.rr,e3.xr,e3.yr);
% fp'rintf(fd,'第4组 %g\t %g\t %g\t %g\t %g\t\n\n\n',e4.evapoints,e4.tr,e4.rr,e4.xr,e4.yr);
% fclose(fd);

display(e1);
display(e2);
display(e3);
display(e4);