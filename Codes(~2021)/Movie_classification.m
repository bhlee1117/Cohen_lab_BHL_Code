%% load 
clear
[fnm,pth]=uigetfile('*.mat','Select all NF_maps','Multiselect','on');
img_pth='E:\BACKUP\대학원\연구실\MY_Projects\In_vivo_imaging\Cell_detection\TXN_detection_deep_learning\Classification\OE9\';
max_file=fullfile(pth,'MAX.tif');
%%
for i=1:length(fnm)
    load([pth fnm{1,i}])
    list(:,i)=NF_map.list(:,4);
    
    sp=split(char(fnm{1,i}),'.');
   sp_tr=split(char(sp{1,1}),'_');
   setDir=[img_pth char(sp_tr{3,1})];
imds = imageDatastore(setDir,'IncludeSubfolders',true,'LabelSource',...
    'foldernames');
Files{i}=imds.Files;
labels{i}=imds.Labels;
labels{i}(labels{i}=='Non_TXN')='No TXN';
   sp=split(imds.Files,'.');
   for j=1:size(sp,1)
   sp_name=split(sp{j,1},'\');
   cell_num(j,i)=str2num(sp_name{13,1});
   end
   
   for j=1:size(list,1)
       row=find(cell_num(:,i)==j);
       cell_file(j,i)=row;
   end
end

%%
aviobj = VideoWriter([pth,'Classification.avi']);
open(aviobj);
maxim=imread(max_file);

for i=1:size(cell_file,1)
    figure(1)
    
    for j=1:length(fnm)
         sp=split(char(fnm{1,j}),'.');
   sp_tr=split(char(sp{1,1}),'_');
   if j>=3 
      subplot(3,length(fnm),j)
   else
   if j==1
    subplot(3,length(fnm),2)
   end
   if j==2
       subplot(3,length(fnm),1)
   end
   end
    im=imread(char(Files{1,j}(cell_file(i,j),1)));
    imagesc(im)
    colormap('gray')
    axis equal
    xlabel(char(labels{1,j}(cell_file(i,j),1)))
    title(char(sp_tr{3,1}))
    xticks([])
    yticks([])
    axis tight
    end
    
    subplot(3,length(fnm),length(fnm)+1:3*length(fnm))
    imagesc(maxim)
    colormap('gray')
    axis equal
    hold all
    plot(NF_map.Cell_position(i,1),NF_map.Cell_position(i,2),'ro','markersize',10,'linewidth',3)
    xticks([])
    yticks([])
    axis tight
    F=figure(1);
writeVideo(aviobj,getframe(F));
if ~mod(i,100)
        close(figure(1))
    end
end
close(figure(1));
 close(aviobj);