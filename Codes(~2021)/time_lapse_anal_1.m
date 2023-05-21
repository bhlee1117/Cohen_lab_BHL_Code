%% Load the image and Lipofuscin tracking data
clear
[fnm,pth]=uigetfile('*.tif','Select the timelapse images','Multiselect','on');
[fnm_lip,pth_lip]=uigetfile('*.tif','Select the timelapse images','Multiselect','on');
[track_fn track_pn]=uigetfile('*.txt','Select the lipofuscin track','Multiselect','on');
channel=2;
zstack=51;
widxz=[50 50];
%%
if ischar(fnm)
    fnm={fnm};
end
tr_org=importdata([track_pn track_fn]);
tr=round(tr_org-tr_org(1,:));

for i=1:size(fnm)
imginf=imfinfo([pth char(fnm{1,i})]);
numstack=size(imginf,1);    
time=numstack/channel/zstack;
end
%%
for i=1:time
    for j=1:zstack
     im3d(:,:,j)=imread([pth fnm{1,1}],2*j-1+(i-1)*channel*zstack);  
     im3d_lip(:,:,j)=imread([pth_lip fnm_lip],2*j-1+(i-1)*channel*zstack);
    end
    imxz(:,:,i)=im3d_lip(round(tr_org(i,2)),:,:);
   
    if i==1
    imwrite(uint16(imxz(:,:,i)),[pth_lip '\Mis\' 'Imxz_' char(fnm_lip)],'Compression','none');
else
imwrite(uint16(imxz(:,:,i)),[pth_lip '\Mis\' 'Imxz_' char(fnm_lip)],'Writemode','Append','Compression','none');
end

end
