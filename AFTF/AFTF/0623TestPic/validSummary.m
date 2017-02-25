function [ out, m, details ] = validSummary( input )
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
m.rthresh=2;
m.xthresh=20;
m.ythresh=20;
m.rbound=20;
m.xbound=50;
m.ybound=50;

details.m=m;

r3d=input(:,3);%rotation error
x3d=input(:,4);%x error
y3d=input(:,5);%y error
rd=input(:,8);
xd=input(:,9);
yd=input(:,10);
fp=input(:,6);%filtered points
cp=input(:,7);%keypoints
pm=input(:,11);%mean dist between each pair points
pv=input(:,12);%variance among all pairs

r3v=0;
x3v=0;
y3v=0;
t=0;
[size3,n]=size(r3d);
out.evapoints=sum(fp)/size3;
f1=0;
f2=0;
f3=0;

rpi=1;
rni=1;
xpi=1;
xni=1;
ypi=1;
yni=1;
tpi=1;
for n=1:size3
    if(abs(r3d(n))<=m.rthresh)
        r3v=r3v+1;
        f1=1;
        details.rp.pm(rpi)=pm(n);
        details.rp.pv(rpi)=pv(n);
        details.rp.fp(rpi)=fp(n);
        details.rp.no(rpi)=n;
        rpi=rpi+1;
    else
        details.rn.pm(rni)=pm(n);
        details.rn.pv(rni)=pv(n);
        details.rn.fp(rni)=fp(n);
        details.rn.no(rni)=n;
        rni=rni+1;
    end
    if(abs(x3d(n))<=m.xthresh)
        x3v=x3v+1;
        f2=1;
        details.xp.pm(xpi)=pm(n);
        details.xp.pv(xpi)=pv(n);
        details.xp.fp(xpi)=fp(n);
        details.xp.no(xpi)=n;
        xpi=xpi+1;
    else
        details.xn.pm(xni)=pm(n);
        details.xn.pv(xni)=pv(n);
        details.xn.fp(xni)=fp(n);
        details.xn.no(xni)=n;
        xni=xni+1;
    end
    if(abs(y3d(n))<=m.ythresh)
        y3v=y3v+1;
        f3=1;
        details.yp.pm(ypi)=pm(n);
        details.yp.pv(ypi)=pv(n);
        details.yp.fp(ypi)=fp(n);
        details.yp.no(ypi)=n;
        ypi=ypi+1;
    else
        details.yn.pm(yni)=pm(n);
        details.yn.pv(yni)=pv(n);
        details.yn.fp(yni)=fp(n);
        details.yn.no(yni)=n;
        yni=yni+1;
    end
    if((f1+f2+f3)==3)
        t=t+1;
        f1=0;
        f2=0;
        f3=0;
        details.tp.pm(tpi)=pm(n);
        details.tp.pv(tpi)=pv(n);
        details.tp.fp(tpi)=fp(n);
        details.tp.no(tpi)=n;
        tpi=tpi+1;
    end
    if(abs(xd(n))>m.xbound || abs(yd(n))>m.ybound || abs(rd(n))>m.rbound)
        t=t+1;
    end
end
out.tr=t/size3;
out.rr=r3v/size3;
out.xr=x3v/size3;
out.yr=y3v/size3;
end

