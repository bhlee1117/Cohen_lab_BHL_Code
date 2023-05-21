%  Image enhancement code
%
%     Website : http://neurobiophysics.snu.ac.kr/
%
% This code is written in purpose of classifying Arc-TXN neurons manually.
% Written by Byung Hun Lee, Deptartment of Physics and Astronomy, Seoul National University, Aug 7th 2020.
% Modification Sep 7 th 2020, Remove filtered montage.
%% Classify the cell in the first and second frame
clear
Class_folder=uigetdir;
[fnm_cell,pth_cell]=uigetfile('.mat','Select the matlab file that contains cell_list and aligned images');
Cell_list_file=fullfile(pth_cell,fnm_cell);
tic;load(Cell_list_file)
toc;
disp('Images and cell list are loaded')
Cell_size=[71 71 71];
Window_size=300;
image_write=false;
       %%
    sp=split(fnm_cell,'.');
    spp=split(sp{1,1},'_');
    sp_name=spp{1,1};
for i=1:size(im_aligned,2)
    mkdir([Class_folder '\' sp_name '\' num2str(i)])
    mkdir([Class_folder '\' sp_name '\' num2str(i) '\TXN'])
    mkdir([Class_folder '\' sp_name '\' num2str(i) '\Non_TXN'])
    mkdir([Class_folder '\' sp_name '\' num2str(i) '\NA'])
end
dir_ts= dir([Class_folder '\' sp_name '\' num2str(1) '\TXN']);
dir_no= dir([Class_folder '\' sp_name '\' num2str(1) '\Non_TXN']);
dir_NA=dir([Class_folder '\' sp_name '\' num2str(1) '\NA']);

dir_pre_cells=[dir_no(3:size(dir_no,1)); dir_ts(3:size(dir_ts,1)) ; dir_NA(3:size(dir_NA,1))];
fnd_result_file=dir(Class_folder);
if ~isempty(cell2mat(strfind({fnd_result_file.name},['NF_map.mat'])))
load([Class_folder,'\NF_map.mat'])
else
    NF_map.list=[];
    ss=randn(1,size(im_aligned,2));
    [ss rewnsa]=sort(ss,'descend');
end
if ~size(dir_pre_cells,1)==0  % Find the date in the categorized folder
    start_here=size( dir_pre_cells,1)+1;
else
    start_here=1;
end
%% Make expanded images
    clear im_exp
    NF_map.Cell_position=cell_list;
    
for i=1:size(im_aligned,2)
    im_exp{rewnsa(1,i)}=zeros(size(im_aligned{i},1)+2*floor(Cell_size(1,1)/2),size(im_aligned{i},2)+2*floor(Cell_size(1,2)/2),size(im_aligned{i},3)+2*floor(Cell_size(1,3)/2));
    im_exp{rewnsa(1,i)}(floor(Cell_size(1,2)/2)+1:floor(Cell_size(1,2)/2)+size(im_aligned{i},1),floor(Cell_size(1,1)/2)+1:floor(Cell_size(1,1)/2)+size(im_aligned{i},2),...
           floor(Cell_size(1,3)/2)+1:floor(Cell_size(1,3)/2)+size(im_aligned{i},3))=im_aligned{i};
    end

    %%
 
    NF_map.Cell_position=cell_list;
n=5;
    for k=1:size(Cell_size,2)
        Cells.list(:,k)=cell_list(:,k)+floor(Cell_size(1,k)/2);
    end
    j=start_here;

    handles.fig = figure(1);
set(gcf, 'Position',  [10, 10, Window_size*size(im_aligned,2)+200, 1150])
cl=find(cell_list(:,3)>12 & cell_list(:,3)<size(im_aligned{1},3)-12);
cl_diff=find(cell_list(:,3)<12 | cell_list(:,3)>size(im_aligned{1},3)-12);
    while j<=size(cl,1)    % cell
    clear crop_im max_crop_im mont filt_mont;
clf('reset')

    for i=1:size(im_aligned,2) % session
crop_im{i}=im_exp{i}(round(Cells.list(cl(j),2))-floor(Cell_size(1,2)/2):round(Cells.list(cl(j),2))+floor(Cell_size(1,2)/2),...
                     round(Cells.list(cl(j),1))-floor(Cell_size(1,1)/2):round(Cells.list(cl(j),1))+floor(Cell_size(1,1)/2),...
                     round(Cells.list(cl(j),3))-floor(Cell_size(1,3)/2):round(Cells.list(cl(j),3))+floor(Cell_size(1,3)/2));
crop_im{i}(crop_im{i}<0)=0;
max_crop_im{i}=max(crop_im{i},[],3);

for jj=1:size(crop_im{i},3) % Filtering image
Bgimfilt=imgaussfilt(crop_im{i}(:,:,jj), 6);
subim=abs(crop_im{i}(:,:,jj)-Bgimfilt);
Filt_im{i}(:,:,jj)=imgaussfilt(subim,1.5);
end

for jj=1:size(crop_im{i},3) % Montage image making
              if mod(jj,n)==0
mont{i}(1+(floor(jj/n)-1)*Cell_size(1,1):(floor(jj/n))*Cell_size(1,1),1+(n-1)*Cell_size(1,2):(n)*Cell_size(1,2))=crop_im{i}(:,:,jj);
% filt_mont{i}(1+(floor(jj/n)-1)*Cell_size(1,1):(floor(jj/n))*Cell_size(1,1),1+(n-1)*Cell_size(1,2):(n)*Cell_size(1,2))=Filt_im{i}(:,:,jj);
           else
mont{i}(1+floor(jj/n)*Cell_size(1,1):(floor(jj/n)+1)*Cell_size(1,1),1+(mod(jj,n)-1)*Cell_size(1,2):(mod(jj,n))*Cell_size(1,2))=crop_im{i}(:,:,jj);
% filt_mont{i}(1+floor(jj/n)*Cell_size(1,1):(floor(jj/n)+1)*Cell_size(1,1),1+(mod(jj,n)-1)*Cell_size(1,2):(mod(jj,n))*Cell_size(1,2))=Filt_im{i}(:,:,jj);
           end
end

%handles.axes1 = axes('Units','pixels','Position',[(Window_size+10)*(i-1)+50 50 Window_size 350]);   

% imshow(filt_mont{i},[min(min(filt_mont{i})) 1.3*max(max(filt_mont{i}))])
% colormap('gray')
% axis equal tight off
% hold on
handles.axes2 = axes('Units','pixels','Position',[(Window_size+10)*(i-1)+50 900 120 120]);
imagesc(max_crop_im{i})
axis equal tight off
colormap('gray')

handles.axes3 = axes('Units','pixels','Position',[(Window_size+10)*(i-1)+50 50 Window_size 850]);
imshow(mont{i},[min(min(mont{i})) 1.1*max(max(mont{i}))])
colormap('gray')
axis equal tight off
handles.axes4 = axes('Units','pixels','Position',[(Window_size+10)*(i-1)+170 900 120 120]);
imagesc(max(Filt_im{i},[],3))
axis equal tight off
colormap('gray')
hold all
    end

    i=1;
    clear sw
while i<=size(im_aligned,2) % session
try
sw(1,i)=input(['( ' num2str(cl(j)) '/' num2str(size(Cells.list,1)) ' )th cell of ' num2str(i) '/' num2str(size(im_aligned,2)) ' th image. Is the cell has TXN site?\n']);
catch
save([Class_folder,'\NF_map.mat'],'NF_map','rewnsa');    
end
if sw(1,i)==4 || sw(1,i)==6
    break
end
if sw(1,i)==5
    i=i-1;
else
    i=i+1;
end
end

if ~isempty(find(sw==6))  % GO back to previous cell
        j=j-1;    
        sw_del=NF_map.list(cl(j),4:end);
for jj=1:size(sw_del,2)

    if sw_del(1,jj)==1
        delete([Class_folder '\' sp_name '\' num2str(jj) '\Non_TXN\',num2str(cl(j)),'.tif'])
    end
    
    if sw_del(1,jj)==2
        delete([Class_folder '\' sp_name '\' num2str(jj) '\TXN\',num2str(cl(j)),'.tif'])
    end
    
    if sw_del(1,jj)==3
        delete([Class_folder '\' sp_name '\' num2str(jj) '\NA\',num2str(cl(j)),'.tif'])
    end
    
end

else  if ~isempty(find(sw==4))  % Break
    close(figure(1))
    close(figure(2))
    break
end

for k=1:size(im_aligned,2) %Save the tif images
if sw(1,k)==1 
    imwrite(uint16(max(crop_im{k}(:,:,jj),[],3)),[Class_folder '\' sp_name '\' num2str(k) '\Non_TXN\',num2str(cl(j)),'.tif'],'Compression','none');
end
if sw(1,k)==2
   imwrite(uint16(max(crop_im{k}(:,:,jj),[],3)),[Class_folder '\' sp_name '\' num2str(k) '\TXN\',num2str(cl(j)),'.tif'],'Compression','none');
end
if sw(1,k)==3
   imwrite(uint16(max(crop_im{k}(:,:,jj),[],3)),[Class_folder '\' sp_name '\' num2str(k) '\NA\',num2str(cl(j)),'.tif'],'Compression','none');
end
end

NF_map.list(cl(j),:)=[NF_map.Cell_position(cl(j),1) NF_map.Cell_position(cl(j),2) Cells.list(cl(j),3)-floor(Cell_size(1,3)/2) sw];
j=j+1;
end
    end

    
    
  close all

  save([Class_folder,'\NF_map.mat'],'NF_map','rewnsa');

%% Need to fix *******2020/08/13
    
    for j=cl_diff'
        for k=1:size(im_aligned,2) %Save the tif images
        imwrite(uint16(zeros(Cell_size(1,1),Cell_size(1,2))),[Class_folder '\' sp_name '\' num2str(k) '\NA\',num2str(j),'.tif'],'Compression','none');
        end
        NF_map.list(j,:)=[NF_map.Cell_position(j,1) NF_map.Cell_position(j,2) Cells.list(j,3)-floor(Cell_size(1,3)/2) zeros(1,size(im_aligned,2))+3];
    end
    
    dir_ts= dir([Class_folder '\' sp_name '\' num2str(1) '\TXN']);
dir_no= dir([Class_folder '\' sp_name '\' num2str(1) '\Non_TXN']);
dir_NA=dir([Class_folder '\' sp_name '\' num2str(1) '\NA']);

dir_pre_cells=[dir_no(3:size(dir_no,1)); dir_ts(3:size(dir_ts,1)) ; dir_NA(3:size(dir_NA,1))];
%%
    
    for f=1:size(im_aligned,2)
if 1
pth=[Class_folder '\' sp_name '\' num2str(f)];
imds = imageDatastore(pth,'IncludeSubfolders',true,'LabelSource',...
    'foldernames');
sw={'Non_TXN','TXN','NA'};

tmp=split(imds.Files,'.');
for ii=1:size(tmp,1)
spp_name{ii,1}=split(tmp{ii,1},'\');
end
for ii=1:size(spp_name{1,1},1)
    if ~isempty(findstr(char(spp_name{1,1}(ii,1)),'NA')) || ~isempty(findstr(char(spp_name{1,1}(ii,1)),'Non_TXN')) || ~isempty(findstr(char(spp_name{1,1}(ii,1)),'TXN'))
        ref=ii;
    end
end
clear tmp

NF=zeros(size(imds.Files,1),1); 
for i=1:size(imds.Files,1)
    for jj=1:size(sw,2)
       if ~isempty(strfind(sw{1,jj},char(imds.Labels(i,1)))) && ~isempty(strfind(char(imds.Labels(i,1)),sw{1,jj}))
           g=str2num(cell2mat(spp_name{i,1}(ref+1,1)));
           if NF(g,1)==0
    NF(g,1)=jj;
           else
               error(['There are two' num2str(g) 'cell in different folder !!!!!!'])
           end
           end
    end   
end
NF_map.list(:,f+3)=NF;
NF_map.list(:,1:3)=NF_map.Cell_position;
else
    error(['Classification is not done yet'])
end

    end
NF_map.mousename=char(sp_name);
save([Class_folder,'\NF_map.mat'],'NF_map','rewnsa');