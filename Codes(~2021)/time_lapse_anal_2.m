clear
[fnm,pth]=uigetfile('*.tif','Select the timelapse images','Multiselect','on');
[track_fn track_pn]=uigetfile('*.txt','Select the lipofuscin track','Multiselect','on');
[track_xz track_pxz]=uigetfile('*.txt','Select the lipofuscin XZ track','Multiselect','on');
channel=2;
zstack=51;
widxz=[50 50];
%%
if ischar(fnm)
    fnm={fnm};
end
tr_org=importdata([track_pn track_fn]);
trxz=importdata([track_pxz track_xz]);
tr_org(:,3)=trxz(:,1);
tr=round(tr_org-tr_org(1,:));

for i=1:size(fnm)
imginf=imfinfo([pth char(fnm{1,i})]);
numstack=size(imginf,1);    
time=numstack/channel/zstack;
end
%%

sub_reg=zeros(imginf(1).Height+abs(min(tr(:,2)))+abs(max(tr(:,2))),imginf(1).Width+abs(min(tr(:,1)))+abs(max(tr(:,1))),zstack+abs(min(tr(:,3)))+abs(max(tr(:,3))));
mask=zeros(imginf(1).Height+abs(min(tr(:,2)))+abs(max(tr(:,2))),imginf(1).Width+abs(min(tr(:,1)))+abs(max(tr(:,1))),zstack+abs(min(tr(:,3)))+abs(max(tr(:,3))))+1;
im=zeros(imginf(1).Height+abs(min(tr(:,2)))+abs(max(tr(:,2))),imginf(1).Width+abs(min(tr(:,1)))+abs(max(tr(:,1))),2);
%%
for i=1:time
maskk=zeros(imginf(1).Height+abs(min(tr(:,2)))+abs(max(tr(:,2))),imginf(1).Width+abs(min(tr(:,1)))+abs(max(tr(:,1))),zstack+abs(min(tr(:,3)))+abs(max(tr(:,3))));
for j=1:zstack
   im=zeros(imginf(1).Height+abs(min(tr(:,2)))+abs(max(tr(:,2))),imginf(1).Width+abs(min(tr(:,1)))+abs(max(tr(:,1))),2);
im(abs(max(tr(:,2)))-tr(i,2)+1:abs(max(tr(:,2)))-tr(i,2)+imginf(1).Height,abs(max(tr(:,1)))-tr(i,1)+1:abs(max(tr(:,1)))-tr(i,1)+imginf(1).Width,1)=imread([pth fnm{1,1}],2*j-1+(i-1)*channel*zstack);
im(abs(max(tr(:,2)))-tr(i,2)+1:abs(max(tr(:,2)))-tr(i,2)+imginf(1).Height,abs(max(tr(:,1)))-tr(i,1)+1:abs(max(tr(:,1)))-tr(i,1)+imginf(1).Width,2)=imread([pth fnm{1,1}],2*j+(i-1)*channel*zstack);
sub_reg(:,:,abs(max(tr(:,3)))-tr(i,3)+j,i)=subtraction_image_ftn(im);

    maskk(:,:,abs(max(tr(:,3)))-tr(i,3)+j)=im(:,:,1)~=0;
end
mask=mask&maskk;
end

%%
[row col]=find(sum(mask,3)>0);
cro_pos=[min(row) min(col) max(row)-min(row) max(col)-min(col)];

for i=1:time
for j=find(sum(sum(mask))>0)'
    
    imwrite(uint16(imcrop(sub_reg(:,:,j,i),cro_pos)),[pth 'subreg_' char(fnm{1,1})],'Writemode','Append','Compression','none');
%imwrite(uint16(mask(:,:,j)),[pth 'mask_' char(fnm{1,1})],'Writemode','Append','Compression','none');
end
end