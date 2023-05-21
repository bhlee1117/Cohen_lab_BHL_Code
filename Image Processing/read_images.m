% 8/1/08 Malcolm Campbell
% Reads in images, reduces to intensity vs time, substracts dark count

numd=50;
numr=500;
numm=200;
nump=200;
numbwlcvr=500;
numbnolcvr=200;
numbwdcm=500;

if 1
dim=read_tifs2('Dark',0,numd-1);
di=squeeze(mean(mean(dim)));
avdi=mean(di);
clear dim;clear di;
end;

% Don't need numd anymore
clear numd;

if 1
ri1=zeros(numr,1);
for j = 0:numr-1
    temp=read_tifs2('RefImage1',j,j);
    ri1(j+1)=mean(mean(temp))-avdi;
    clear temp;
end;
end;

if 1
ri2=zeros(numr,1);
for j = 0:numr-1
    temp=read_tifs2('RefImage2',j,j);
    ri2(j+1)=mean(mean(temp))-avdi;
    clear temp;
end;
end;

% Don't need numr anymore
clear numr;

if 1
mi1=zeros(numm,1);
for j = 0:numm-1
    temp=read_tifs2('m2Image1',j,j);
    mi1(j+1)=mean(mean(temp))-avdi;
    clear temp;
end;
end;

if 1
mi2=zeros(numm,1);
for j = 0:numm-1
    temp=read_tifs2('m2Image2',j,j);
    mi2(j+1)=mean(mean(temp))-avdi;
    clear temp;
end;
end;

% Don't need numm anymore
clear numm;

if 1
pi1=zeros(nump,1);
for j = 0:nump-1
    temp=read_tifs2('p2Image1',j,j);
    pi1(j+1)=mean(mean(temp))-avdi;
    clear temp;
end;
end;

if 1
pi2=zeros(nump,1);
for j = 0:nump-1
    temp=read_tifs2('p2Image2',j,j);
    pi2(j+1)=mean(mean(temp))-avdi;
    clear temp;
end;
end;

% Don't need nump anymore
clear nump;

if 0
bwlcvri1=zeros(numbwlcvr,1);
for j = 0:numbwlcvr-1
    temp=read_tifs2('blankwLCVR1',j,j);
    bwlcvri1(j+1)=mean(mean(temp))-avdi;
    clear temp;
end;
end;

if 0
bwlcvri2=zeros(numbwlcvr,1);
for j = 0:numbwlcvr-1
    temp=read_tifs2('blankwLCVR2',j,j);
    bwlcvri2(j+1)=mean(mean(temp))-avdi;
    clear temp;
end;
end;

% Don't need numbwlcvr anymore
clear numbwlcvr

if 0
bnolcvri1=zeros(numbnolcvr,1);
for j = 0:numbnolcvr-1
    temp=read_tifs2('blanknoLCVR1',j,j);
    bnolcvri1(j+1)=mean(mean(temp))-avdi;
    clear temp;
end;
end;

if 0
bnolcvri2=zeros(numbnolcvr,1);
for j = 0:numbnolcvr-1
    temp=read_tifs2('blanknoLCVR2',j,j);
    bnolcvri2(j+1)=mean(mean(temp))-avdi;
    clear temp;
end;
end;

% Don't need numbnolcvr anymore
clear numbnolcvr;

if 0
bwdcmi1=zeros(numbwdcm,1);
for j = 0:numbwdcm-1
    temp=read_tifs2('blankwDCM1',j,j);
    bwdcmi1(j+1)=mean(mean(temp))-avdi;
    clear temp;
end;
end;

if 0
bwdcmi2=zeros(numbwdcm,1);
for j = 0:numbwdcm-1
    temp=read_tifs2('blankwDCM2',j,j);
    bwdcmi2(j+1)=mean(mean(temp))-avdi;
    clear temp;
end;
end;

% Don't need numbwdcm anymore
clear numbwdcm;

clear j;