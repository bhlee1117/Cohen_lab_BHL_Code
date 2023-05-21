%  Classification code for RNA scope image.
%
%     Website : http://neurobiophysics.snu.ac.kr/
%
% INPUTS 

% RNA scope image

% OUTPUTS
% 1. Cropped cell image
% 2. TXN list
%
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy, Seoul National University, 11/20/2018.
%% Classify the cell in the first and second frame
clear
[fnm pth] = uigetfile('.mat','Select the mat file which contain the variable RNAscope','Multiselect','on'); %
[fnmptl pthptl] = uigetfile('.mat','Select the trackNtrace file which contain the information about c-fos particles','Multiselect','on'); %
save_dir='H:\Image_data\KIST_RNAscope\20191101_SNU-RNAscope_Arc-Fos-Egr1\Analysis\';

Class_folder='H:\Image_data\KIST_RNAscope\20191101_SNU-RNAscope_Arc-Fos-Egr1\Analysis\Classification\';
Cell_size=[151 151];
expander=301;
ch=4;
%% 

clear start_here
    sp=split(fnm,'.');
    load([pth fnm])
    load([pthptl fnmptl])
   for i=1:ch-1
    mkdir([Class_folder sp{1,1} '\' num2str(i) '\TXN'])
    mkdir([Class_folder sp{1,1} '\' num2str(i) '\Non_TXN'])
    mkdir([Class_folder sp{1,1} '\' num2str(i) '\NA'])
   end
   
  


%%
    clear im im_exp
    
    im_exp=zeros(size(RNAscope.max,1)+2*floor(Cell_size(1,1)/2),size(RNAscope.max,2)+2*floor(Cell_size(1,2)/2),4);
    im_exp(floor(Cell_size(1,2)/2)+1:floor(Cell_size(1,2)/2)+size(RNAscope.max,1),floor(Cell_size(1,1)/2)+1:floor(Cell_size(1,1)/2)+size(RNAscope.max,2),1:4)=RNAscope.max;
    
    im_tre_exp=zeros(size(RNAscope.max,1)+2*floor(expander/2),size(RNAscope.max,2)+2*floor(expander/2),4);
    im_tre_exp(floor(expander/2)+1:floor(expander/2)+size(RNAscope.max,1),floor(expander/2)+1:floor(expander/2)+size(RNAscope.max,2),1:4)=RNAscope.max;
    RNAscope.result=[];
  
    %%
   
    start_here=size(RNAscope.result,1)+1;
    color=[0.7 .7 .7;0 .5 0;1.1 0 0;0 0.4 0.4];
    sigma_thres=2.2;
    amp_thres=60;
          sp=split(fnm,'.');
          clear Cells.list
           for k=1:size(RNAscope.list,2) 
        Cells.list(:,k)=RNAscope.list(:,k)+floor(Cell_size(1,k)/2);
        Cells.list_exp(:,k)=RNAscope.list(:,k)+floor(expander/2);
           end
           j=start_here;
while j<=size(RNAscope.list,1) % cells 

    clear crop_im max_crop_im mont filt_mont RGB_im RGB_im_exp;

    crop_im=im_exp(round(Cells.list(j,2))-floor(Cell_size(1,2)/2):round(Cells.list(j,2))+floor(Cell_size(1,2)/2),... % make crop images
               round(Cells.list(j,1))-floor(Cell_size(1,1)/2):round(Cells.list(j,1))+floor(Cell_size(1,1)/2),1:4);
           
    crop_exp=im_tre_exp(round(Cells.list_exp(j,2))-floor(expander/2):round(Cells.list_exp(j,2))+floor(expander/2),...
               round(Cells.list_exp(j,1))-floor(expander/2):round(Cells.list_exp(j,1))+floor(expander/2),1:4);
    ptls_fd=find((refinementData{1}(:,2)>round(RNAscope.list(j,2))-floor(Cell_size(1,2)/2) & refinementData{1}(:,2)<round(RNAscope.list(j,2))+floor(Cell_size(1,2)/2)...
               & refinementData{1}(:,1)>round(RNAscope.list(j,1))-floor(Cell_size(1,1)/2) & refinementData{1}(:,1)<round(RNAscope.list(j,1))+floor(Cell_size(1,1)/2)...
               & refinementData{1}(:,6).*refinementData{1}(:,7)>sigma_thres & refinementData{1}(:,4)>amp_thres)==1);
    ptls=refinementData{1}(ptls_fd,:);
    RGB_im=zeros(size(crop_im,1),size(crop_im,2),3);
    RGB_im_exp=zeros(size(crop_exp,1),size(crop_exp,2),3);
    for col=1:size(color,1)
        for rgb=1:size(color,2)
    RGB_im(:,:,rgb)=RGB_im(:,:,rgb)+crop_im(:,:,col)/max(max(crop_im(:,:,col)))*color(col,rgb);
    RGB_im_exp(:,:,rgb)=RGB_im_exp(:,:,rgb)+crop_exp(:,:,col)/max(max(crop_exp(:,:,col)))*color(col,rgb);
        end
    end
    
% RGB_im=flipud(RGB_im);
% RGB_im_exp=flipud(RGB_im_exp);% Showing popup window
    figure(1)
subplot(3,2,1)
imagesc(RGB_im)
axis equal tight
subplot(3,2,2)
imagesc(RGB_im_exp)
colormap('gray')
axis equal tight
hold all
plot(expander/2,expander/2,'go','markersize',70)
subplot(3,2,3) % Arc
imagesc(crop_im(:,:,2))
axis equal tight
subplot(3,2,4) % c-fos
imagesc(crop_im(:,:,3))
hold all
plot(ptls(:,1)-RNAscope.list(j,1)+floor(Cell_size(1,1)/2),ptls(:,2)-RNAscope.list(j,2)+floor(Cell_size(1,2)/2),'ro','markersize',12)
axis equal tight
subplot(3,2,5)
imagesc(crop_im(:,:,4))
axis equal tight
i=1;
while i<=ch-1  %switch cells

sw(1,i)=input(['( ' num2str(j) '/' num2str(size(Cells.list,1)) ' )th cell of ' num2str(i) '/' num2str(ch-1) ' th image. Is the cell has TXN site?\n']);

if sw(1,i)==4 || sw(1,i)==6
    break
end
if sw(1,i)==5   % Back to previous channel
    i=i-1;
else
    i=i+1;
end
end


if ~isempty(find(sw==6))
        j=j-1;    
        sw_del=RNAscope.result(j,4:end);
for jj=1:size(sw_del,2) %channel

    if sw_del(1,jj)==1
        delete([Class_folder sp{1,1} '\' num2str(jj) '\Non_TXN\',num2str(j),'.tif'])
    end
    
    if sw_del(1,jj)==2
        delete([Class_folder sp{1,1} '\' num2str(jj) '\TXN\',num2str(j),'.tif'])
    end
    
    if sw_del(1,jj)==3
        delete([Class_folder sp{1,1} '\' num2str(jj) '\NA\',num2str(j),'.tif'])
    end
    
end



else  if ~isempty(find(sw==4))
    close(figure(1))
    close(figure(2))
    break
end

for k=1:ch-1 %channel

if sw(1,k)==1
    imwrite(uint16(crop_im(:,:,k+1)),[Class_folder sp{1,1} '\' num2str(k) '\Non_TXN\',num2str(j),'.tif'],'Compression','none');
end

if sw(1,k)==2
    imwrite(uint16(crop_im(:,:,k+1)),[Class_folder sp{1,1} '\' num2str(k) '\TXN\',num2str(j),'.tif'],'Compression','none');
end

if sw(1,k)==3
    imwrite(uint16(crop_im(:,:,k+1)),[Class_folder sp{1,1} '\' num2str(k) '\NA\',num2str(j),'.tif'],'Compression','none');
end
end
RNAscope.result(j,:)=[RNAscope.list(j,1) RNAscope.list(j,2) sw];
j=j+1;
end
end
  close(figure(1))
    close(figure(2))
save([Class_folder,'NF_map_',sp{1,1},'.mat'],'RNAscope');


