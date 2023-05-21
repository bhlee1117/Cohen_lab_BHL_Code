%% Load registered image
clear

[ fnm pth]=uigetfile('*.tif','Select the TXN data','Multiselect','on');
[fnm_ca pth_ca ]=uigetfile('*.tif','Select the Calcium image');
fnm_ca={fnm_ca};
x_pixel=0.18; %um
z_pixel=0.25; %um
ratio=x_pixel/z_pixel;
%% imaging loading
    clear imginf im_TXN im_ca
imginf=imfinfo([pth fnm{1,1}]);
numstack=size(imginf,1);
for i=1:numstack
    for j=1:length(fnm)
    im_TXN{j}(:,:,i)=imread([pth fnm{1,j}],i);
    end
end

imginf=imfinfo([pth_ca fnm_ca{1,1}]);
numstac=size(imginf,1);

for i=1:numstac
    im_ca(:,:,i)=imread([pth_ca fnm_ca{1,1}],i);
end

resize_ratio=size(im_ca,1)/size(im_TXN{1},1);
%% Cell position matching

handles.fig = figure(1);
handles.axes1 = axes('Units','pixels','Position',[50 350 450 450]);
handles.Image = max(max(mean(im_ca,3)))-mean(im_ca,3);
imagesc(handles.Image);
colormap('gray')
axis equal
axis off
for n=1:length(fnm)
handles.axes2= axes('Units','pixels','Position',[500 350 450 450]);
handles.Image2 = im_TXN{n}(:,:,1);
handles.slider = uicontrol('Style','slider','Position',[50 50 500 20],'Min',1,'Max',size(im_TXN{n},3),'Value',3);
handles.Listener = addlistener(handles.slider,'Value','PostSet',@(s,e) im_recall(handles,im_TXN{n}));
imagesc(handles.Image2);
axis equal
axis off
guidata(handles.fig);

zmat(1,n)=input('slider value\n');
imca_ave=mean(imresize(im_ca,1/resize_ratio),3);

ImageToWarp = im_TXN{n}(:,:,zmat);
BaseImage = imca_ave;

tempWarp = double(ImageToWarp);
tempBase = double(BaseImage);

tempWarp = tempWarp./max(max(tempWarp))*3;
tempBase = tempBase./max(max(tempBase));

h = cpselect(tempWarp, tempBase);
uiwait(msgbox('Click OK after closing the CPSELECT window.','Waiting...'))

movingPoints = cpcorr(movingPoints, fixedPoints,tempWarp,tempBase);
tform2 = maketform('projective',movingPoints, fixedPoints);

for z=1:numstack
    ImageToWarp=im_TXN{n}(:,:,z);
[transVis xdata ydata] = imtransform(ImageToWarp, tform2);
transVisWithOffset = imtransform(ImageToWarp, tform2, 'XData', [1 (size(ImageToWarp,2)+tform2.tdata.T(3,1))],'YData', [1 (size(ImageToWarp,1)+ tform2.tdata.T(3,2))]);
transVisBuffrd = uint16(zeros(size(BaseImage)));
transVisBuffrd(1:size(transVisWithOffset,1),1:size(transVisWithOffset,2)) = transVisWithOffset;
Transformation_Matrix=tform2.tdata.T(1:2,1:2);
Translation_Vector=tform2.tdata.T(3,1:2)';
transVisBuffrd = transVisBuffrd((1:size(ImageToWarp,1)), (1:size(ImageToWarp,2)));
im_TXN_reg{n}(:,:,z)=transVisBuffrd;
end

mat_im=zeros(size(imca_ave,1),size(imca_ave,2),3);
close(figure)
end
mat_im(:,:,1)=imca_ave/max(max(imca_ave));
log_TXN_reg=max(im_TXN_reg{1},[],3)>0;
for i=2:length(fnm)
log_TXN_reg=log_TXN_reg&max(im_TXN_reg{i},[],3)>0;
end
mat_im(:,:,2)=log_TXN_reg;
%%
imshow(mat_im)
h = imrect(gca,[10 10 100 100]);
position=wait(h);
crop_pos=round(position);
%%

clear subt_crop res_im reslice im_TXN_reg_crop im_ca_crop
for n=1:length(fnm)
for i=1:81
     im_TXN_reg_crop{n}(:,:,i)=imcrop(im_TXN_reg{n}(:,:,i),crop_pos);
end
end

for i=1:size(im_ca,3)
    im_ca_crop(:,:,i)=imcrop(imresize(im_ca(:,:,i),1/resize_ratio),crop_pos);
end
for n=1:length(fnm)
for i=1:size(im_TXN_reg_crop{n},3)
    sp=split(fnm{1,n},'.');
    spp=split(sp{1,1},'_');
    imwrite(im_TXN_reg_crop{n}(:,:,i),[pth char(spp{end,1}) 'TXN_reg_crop.tif'],'tif','Writemode','Append','Compression','none')
end
end
for i=1:size(im_ca,3)
    imwrite(im_ca_crop(:,:,i),[pth 'Ca_reg_crop.tif'],'tif','Writemode','Append','Compression','none')
end

%%

cell_list=G_cell_find(im_TXN_reg{2},ratio);

%% Slider to set the range
cell_plane=cell_list(find(round(cell_list(:,3))==zmat(1,end)),1:3);
handles.fig = figure;
handles.axes1 = axes('Units','pixels','Position',[55 5 450 450]);
handles.Image = mean(im_ca_crop,3);
handles.cell=cell_plane;
handles.resize_ratio=1;
%rsz=size(max_res_im,1)/size(im_ca_crop,1);
imagesc(handles.Image);
colormap('gray')
axis equal
axis off
hold all

plot(cell_plane(:,1)/1,cell_plane(:,2)/1,'ro','markersize',13)

handles.slider = uicontrol('Style','slider','Position',[20 20 500 20],'Min',1,'Max',80,'Value',1);
handles.Listener = addlistener(handles.slider,'Value','PostSet',@(s,e) pos_recall(handles,cell_list,zmat));
guidata(handles.fig);

%% Crop the ROI
close(figure(1))
range=input('slider value\n');
cell_plane=cell_list(find(round(cell_list(:,3))>=zmat-range & round(cell_list(:,3))<=zmat+range),1:3);
Cell_size=[61 61 61];
%ca_Cell_size=[121 121];
for n=1:length(fnm)
    sp_fnm=split(fnm{1,n},'.');
    sp_fnm=split(sp_fnm{1,1},'_');
    mkdir([pth  '\Classification\' char(sp_fnm{end,1})])
    mkdir([pth  '\Classification\' char(sp_fnm{end,1}) '\TXN'])
    mkdir([pth '\Classification\' char(sp_fnm{end,1}) '\Non_TXN'])
    mkdir([pth '\Classification\' char(sp_fnm{end,1}) '\NA'])
end
%%
for n=1:length(fnm)
    im_exp{n}=zeros(size(im_TXN_reg_crop{n},1)+2*floor(Cell_size(1,1)/2),size(im_TXN_reg_crop{n},2)+2*floor(Cell_size(1,2)/2),size(im_TXN_reg_crop{n},3)+2*floor(Cell_size(1,3)/2));
    im_exp{n}(floor(Cell_size(1,2)/2)+1:floor(Cell_size(1,2)/2)+size(im_TXN_reg_crop{n},1),floor(Cell_size(1,1)/2)+1:floor(Cell_size(1,1)/2)+size(im_TXN_reg_crop{n},2),...
           floor(Cell_size(1,3)/2)+1:floor(Cell_size(1,3)/2)+size(im_TXN_reg_crop{n},3))=im_TXN_reg_crop{n};
end
    %%
    
    n=12;
    if isfile([pth,'NF_map_.mat'])
        load([pth,'NF_map_.mat'])
        start_here=size(NF_map.list,1)+1;
        
    else
        start_here=1;
    end
    for k=1:size(cell_plane,2)
        if k~=3
        cell_plane_shift(:,k)=cell_plane(:,k)+floor(Cell_size(1,k)/2);
        else
            cell_plane_shift(:,k)=cell_plane(:,k)+floor(Cell_size(1,k)/2);
        end
    end
    for nn=1:length(fnm)
        
         sp_fnm=split(fnm{1,nn},'.');
    sp_fnm=split(sp_fnm{1,1},'_');
    
for j=start_here:size(cell_plane_shift,1)    
    clear crop_im max_crop_im mont filt_mont;

    crop_im=im_exp{nn}(round(cell_plane_shift(j,2))-floor(Cell_size(1,2)/2):round(cell_plane_shift(j,2))+floor(Cell_size(1,2)/2),...
               round(cell_plane_shift(j,1))-floor(Cell_size(1,1)/2):round(cell_plane_shift(j,1))+floor(Cell_size(1,1)/2),...
               round(cell_plane_shift(j,3))-floor(Cell_size(1,3)/2):round(cell_plane_shift(j,3))+floor(Cell_size(1,3)/2));

max_crop_im=max(crop_im,[],3);

for jj=1:size(crop_im,3)
Bgimfilt=imgaussfilt(crop_im(:,:,jj), 6);
subim=abs(crop_im(:,:,jj)-Bgimfilt);
Filt_im(:,:,jj)=imgaussfilt(subim,1.5);
end

for jj=1:size(crop_im,3)
              if mod(jj,n)==0
mont(1+(floor(jj/n)-1)*Cell_size(1,1):(floor(jj/n))*Cell_size(1,1),1+(n-1)*Cell_size(1,2):(n)*Cell_size(1,2))=crop_im(:,:,jj);
filt_mont(1+(floor(jj/n)-1)*Cell_size(1,1):(floor(jj/n))*Cell_size(1,1),1+(n-1)*Cell_size(1,2):(n)*Cell_size(1,2))=Filt_im(:,:,jj);
           else
mont(1+floor(jj/n)*Cell_size(1,1):(floor(jj/n)+1)*Cell_size(1,1),1+(mod(jj,n)-1)*Cell_size(1,2):(mod(jj,n))*Cell_size(1,2))=crop_im(:,:,jj);
filt_mont(1+floor(jj/n)*Cell_size(1,1):(floor(jj/n)+1)*Cell_size(1,1),1+(mod(jj,n)-1)*Cell_size(1,2):(mod(jj,n))*Cell_size(1,2))=Filt_im(:,:,jj);
           end
end

figure(1)
subplot(2,3,1:2)
imagesc(mont)
% imshow(mont,[0 max(max(mont))])
colormap('gray')
axis equal
hold on
title('Subtracted Image')
subplot(2,3,3)
imagesc(max_crop_im)
axis equal
colormap('gray')

subplot(2,3,4:5)
imagesc(filt_mont)
colormap('gray')
axis equal
title('Filtered Image')
subplot(2,3,6)
imagesc(max(Filt_im,[],3))
axis equal
colormap('gray')


hold all
sw=input(['( ' num2str(j) '/' num2str(size(cell_plane_shift,1)) ' )th cell. Is the cell has TXN site?\n']);
if sw==1
    for jj=1:size(crop_im,3)
        if jj==1
            imwrite(uint16(crop_im(:,:,jj)),[pth  '\Classification\' char(sp_fnm{end,1}) '\Non_TXN\',num2str(j),'.tif'],'Compression','none');
        else
   imwrite(uint16(crop_im(:,:,jj)),[pth  '\Classification\' char(sp_fnm{end,1}) '\Non_TXN\',num2str(j),'.tif'],'Writemode','Append','Compression','none');
        end
    end
end

if sw==2
   for jj=1:size(crop_im,3)
       if jj==1
           imwrite(uint16(crop_im(:,:,jj)),[pth  '\Classification\' char(sp_fnm{end,1}) '\TXN\',num2str(j),'.tif'],'Compression','none');
       else
   imwrite(uint16(crop_im(:,:,jj)),[pth  '\Classification\' char(sp_fnm{end,1}) '\TXN\',num2str(j),'.tif'],'Writemode','Append','Compression','none');
       end
   end
end

if sw==3
   for jj=1:size(crop_im,3)
       if jj==1
   imwrite(uint16(crop_im(:,:,jj)),[pth  '\Classification\' char(sp_fnm{end,1}) '\NA\',num2str(j),'.tif'],'Compression','none');        
       else
   imwrite(uint16(crop_im(:,:,jj)),[pth  '\Classification\' char(sp_fnm{end,1}) '\NA\',num2str(j),'.tif'],'Writemode','Append','Compression','none');
       end
   end
end
if sw==4
    close(figure(1))
    close(figure(2))
    break
end
NF_map.list(j,[1:3 nn+3])=[cell_plane(j,1) cell_plane(j,2) cell_plane(j,3) sw];
if mod(j,50)==0
close(figure(1))
    close(figure(2))
end
end
end
  close(figure(1))
    close(figure(2))
    
NF_map.cell=cell_plane;
save([pth,'NF_map_.mat'],'NF_map');



