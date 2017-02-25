clear all
threshold=1;
% G1=load('D:\YECWorkspace\CODES\AFTF\AFTF\0623TestPic\G1_20160625_2220_Ransac4\_summary');
% G2=load('D:\YECWorkspace\CODES\AFTF\AFTF\0623TestPic\G2_20160625_2221_Ransac4\_summary');
% G3=load('D:\YECWorkspace\CODES\AFTF\AFTF\0623TestPic\G3_20160625_2222_Ransac4\_summary');
% G4=load('D:\YECWorkspace\CODES\AFTF\AFTF\0623TestPic\G4_20160625_2225_Ransac4\_summary');
% save('G1-G4v4.mat','G1','G2','G3','G4');
load('G1-G4v4.mat');
%%
[rows,cols]=size(G1);
clear G FG
g=1;fg=1;
for i=1:rows
    if(abs(G1(i,9))<threshold)
        G(g,:)=G1(i,:);
        g=g+1;
    else
        FG(fg,:)=G1(i,:);
        fg=fg+1;
    end
end
clear G1
G1=G;
%%
[rows,cols]=size(G2);
clear G
g=1;
for i=1:rows
    if(abs(G2(i,9))<threshold)
        G(g,:)=G2(i,:);
        g=g+1;
    else
        FG(fg,:)=G2(i,:);
        fg=fg+1;
    end
end
clear G2
G2=G;
%%
[rows,cols]=size(G3);
clear G
g=1;
for i=1:rows
    if(abs(G3(i,9))<threshold)
        G(g,:)=G3(i,:);
        g=g+1;
    else
        FG(fg,:)=G3(i,:);
        fg=fg+1;
    end
end
clear G3
G3=G;
%%
[rows,cols]=size(G4);
clear G
g=1;
for i=1:rows
    if(abs(G4(i,9))<threshold)
        G(g,:)=G4(i,:);
        g=g+1;
    else
        FG(fg,:)=G4(i,:);
        fg=fg+1;
    end
end
clear G4
G4=G;
clear G g fg rows cols i
G=[G1;G2;G3;G4];
clear G1 G2 G3 G4
%%
rotation=G(:,3)-G(:,2);
x=G(:,5)-G(:,4);
y=G(:,7)-G(:,6);
ree=G(:,9);
sinv=G(:,20);
cosv=G(:,19);
scale=G(:,21);
cf=G(:,13);
vx=G(:,10);
vy=G(:,11);
convin=G(:,12);
ax=G(:,17);
ay=G(:,18);
c=[G(:,1)';G(:,22)';G(:,23)']';
r_vs_cf=[ree';cf']';
ro_vs_vx_vy=[ree';rotation';x';vx';y';vy']';
pmean=G(:,15);
pvar=G(:,16);
xr=G(:,25);
yr=G(:,26);
disvar=G(:,24);
ro_vs_pmpv_cf=[rotation';ree';pmean';pvar';cf';disvar']';
rx=x-xr;
ry=y-yr;
cosvx=cosv.*vx;
cosvax=cosv.*ax;
cosvy=cosv.*vy;
cosvay=cosv.*ay;
sinvx=sinv.*vx;
sinvax=sinv.*ax;
sinvy=sinv.*vy;
sinvay=sinv.*ay;
rxy=[rx';ry']';

%%
frotation=FG(:,3)-FG(:,2);
fx=FG(:,5)-FG(:,4);
fy=FG(:,7)-FG(:,6);
free=FG(:,9);
fsinv=FG(:,20);
fcosv=FG(:,19);
fscale=FG(:,21);
fcf=FG(:,13);
fvx=FG(:,10);
fvy=FG(:,11);
fconvin=FG(:,12);
fax=FG(:,17);
fay=FG(:,18);
fdisvar=FG(:,24);
fc=[FG(:,1)';FG(:,22)';FG(:,23)']';
fr_vs_cf=[free';fcf']';
fro_vs_vx_vy=[free';frotation';fx';fvx';fy';fvy']';
fpmean=FG(:,15);
fpvar=FG(:,16);
fxr=FG(:,25);
fyr=FG(:,26);
fro_vs_pmpv_cf=[frotation';free';fpmean';fpvar';fcf';fdisvar']';