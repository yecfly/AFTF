function [ out, m] = validSummary7( input )
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
m.rthresh=3;
m.xthresh=30;
m.ythresh=30;

r3d=input(:,5);%rotation error
x3d=input(:,6);%x error
y3d=input(:,7);%y error

r3v=0;
x3v=0;
y3v=0;
t=0;
[size3,n]=size(r3d);
f1=0;
f2=0;
f3=0;

for n=1:size3
    if(abs(r3d(n))<=m.rthresh)
        r3v=r3v+1;
        f1=1;
    end
    if(abs(x3d(n))<=m.xthresh)
        x3v=x3v+1;
        f2=1;
    end
    if(abs(y3d(n))<=m.ythresh)
        y3v=y3v+1;
        f3=1;
    end
    if((f1+f2+f3)==3)
        t=t+1;
        f1=0;
        f2=0;
        f3=0;
    end
end
out.tr=t/size3;
out.rr=r3v/size3;
out.xr=x3v/size3;
out.yr=y3v/size3;
end

