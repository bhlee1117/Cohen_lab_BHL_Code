%% Load the require images and information
% Cell acceptance code1
% Input : Raw images, Cell coordinate file, zstack, Transformation matrix
clear

[fnm pth]=uigetfile('.tif','Select the raw images','Multiselect','on');
[fnmtr,pthtr]=uigetfile('.mat','Select the tracking file','Multiselect','on');
[fnm_T pth_T]=uigetfile('.mat','Select the transformation matrix','Multiselect','on');
[fnm_C pth_C]=uigetfile('.mat','Select the list of cell','Multiselect','on');
zstack_rng=[61 81;18 38;26 46;1 21;34 54];
%%

if ischar(fnmtr)
    fnmtr={fnmtr};
end
h=waitbar(0,'Drift calculation');

for i=1:length(fnm)
imginf=imfinfo([pth char(fnm{1,i})]);
numstack=size(imginf,1);  
ss=split(fnm{1,i},'.');
sp_name=split(ss{1,1},'_');
g=1;
clear lip_cred_tr lip_cred_tr_cov
for k=1:length(fnmtr)
    if ~isempty(findstr(fnmtr{1,k},sp_name{3,1})) && ~isempty(findstr(sp_name{3,1},fnmtr{1,k}))
%     if ~isempty(findstr(fnmtr{1,k},sp_name{3,1})) && ~isempty(findstr(sp_name{3,1},fnmtr{1,k})) && ~isempty(findstr(fnmtr{1,k},sp_name{6,1})) && ~isempty(findstr(sp_name{6,1},fnmtr{1,k}))...
%        && ~isempty(findstr(fnmtr{1,k},sp_name{2,1})) && ~isempty(findstr(sp_name{2,1},fnmtr{1,k}))
       lip_tr=importdata([pthtr fnmtr{1,k}]);
              lip_cred_tr{1}=zeros(81,1);
       lip_cred_tr{2}=zeros(81,1);
       for ptl=1:max(lip_tr.trackingData(:,1))
        
           if size(find(lip_tr.trackingData(:,1)==ptl),1)>10
               clear dif row
               lip_cred_tr{1}(lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),2),g)=[NaN ; diff(lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),3))];
               lip_cred_tr{2}(lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),2),g)=[NaN ; diff(lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),4))];
                g=g+1;
           end
       end
       lip_cred_tr_cov=lip_cred_tr;
       lip_cred_tr_cov{1}(lip_cred_tr{1}==0)=NaN;
       lip_cred_tr_cov{2}(lip_cred_tr{2}==0)=NaN;
       
       lip_cred_tr_cov{1}(isnan(lip_cred_tr{1}))=0;
       lip_cred_tr_cov{2}(isnan(lip_cred_tr{2}))=0;
       drift{i,1}(1:numstack/2,:)=round(cumsum(mean(lip_cred_tr_cov{1}(1:numstack/2,:),2,'omitnan'),'omitnan'));
       drift{i,2}(1:numstack/2,:)=round(cumsum(mean(lip_cred_tr_cov{2}(1:numstack/2,:),2,'omitnan'),'omitnan'));
    end
end
waitbar((i)/length(fnm))
end
close(h)

h=waitbar(0,'Registration');

for i=1:length(fnm)
    clear im2 subtraction
imginf=imfinfo([pth char(fnm{1,i})]);
numstack=size(imginf,1);  
im2=zeros(imginf(1).Height+abs(min(drift{i,2}(:,1)))*(1-heaviside(min(drift{i,2}(:,1))))+abs(max(drift{i,2}(:,1))),imginf(1).Width+abs(min(drift{i,1}(:,1)))*(1-heaviside(min(drift{i,1}(:,1))))+abs(max(drift{i,1}(:,1))),numstack/2);
cro_pos=[abs(max(drift{i,1}(:,1)))-drift{i,2}(1,1)+1 abs(max(drift{i,2}(:,1)))-drift{i,1}(1,1)+1 imginf(1).Height-1 imginf(1).Width-1];
for j=1:numstack/2   

im=zeros(imginf(1).Height+abs(min(drift{i,2}(:,1)))*(1-heaviside(min(drift{i,2}(:,1))))+abs(max(drift{i,2}(:,1))),imginf(1).Width+abs(min(drift{i,1}(:,1)))*(1-heaviside(min(drift{i,1}(:,1))))+abs(max(drift{i,1}(:,1))),1);
im(abs(max(drift{i,2}(:,1)))-drift{i,2}(j,1)+1:abs(max(drift{i,2}(:,1)))-drift{i,2}(j,1)+imginf(1).Height,...
   abs(max(drift{i,1}(:,1)))-drift{i,1}(j,1)+1:abs(max(drift{i,1}(:,1)))-drift{i,1}(j,1)+imginf(1).Width,1)=imread([pth char(fnm{1,i})],j);

im2(abs(max(drift{i,2}(:,1)))-drift{i,2}(j,1)+1:abs(max(drift{i,2}(:,1)))-drift{i,2}(j,1)+imginf(1).Height,...
   abs(max(drift{i,1}(:,1)))-drift{i,1}(j,1)+1:abs(max(drift{i,1}(:,1)))-drift{i,1}(j,1)+imginf(1).Width,j)=imread([pth char(fnm{1,i})],j);

reg_im{i}(:,:,j)=uint16(imcrop(im2(:,:,j),cro_pos));
end   
waitbar((i)/length(fnm))
end
close(h)
%%
clear im im2
load([pth_C fnm_C])
g=1;
al_crop_im{1}=[];
for j=zstack_rng(1,1):zstack_rng(1,2)
al_crop_im{1}(:,:,g)=imcrop(reg_im{1}(:,:,j),Cells.crop);
imwrite(uint16(al_crop_im{1}(:,:,g)),[pth 'Li_' char(fnm{1,1})],'Writemode','Append','Compression','none');
g=g+1;
end

for i=2:length(fnm)
    load([pth_T fnm_T{1,i-1}])
g=1;
al_crop_im{i}=[];
for z=zstack_rng(i,1):zstack_rng(i,2)
[transVis xdata ydata] = imtransform(reg_im{i}(:,:,z), tform2);
transVisWithOffset = imtransform(reg_im{i}(:,:,z), tform2, 'XData', [1 (size(reg_im{i}(:,:,z),2)+tform2.tdata.T(3,1))],'YData', [1 (size(reg_im{i}(:,:,z),1)+ tform2.tdata.T(3,2))]);
transVisBuffrd = uint16(zeros(size(reg_im{i}(:,:,z))));
transVisBuffrd(1:size(transVisWithOffset,1),1:size(transVisWithOffset,2)) = transVisWithOffset;
Transformation_Matrix=tform2.tdata.T(1:2,1:2);
Translation_Vector=tform2.tdata.T(3,1:2)';
al_crop_im{i}(:,:,g) = imcrop(transVisBuffrd((1:size(reg_im{i}(:,:,z),1)), (1:size(reg_im{i}(:,:,z),2))),Cells.crop);
imwrite(uint16(al_crop_im{i}(:,:,g)),[pth 'Li_' char(fnm{1,i})],'Writemode','Append','Compression','none');
g=g+1;
end

end
%save([pth 'Lip_image.mat'],'al_crop_im')