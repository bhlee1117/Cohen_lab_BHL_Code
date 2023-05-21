%% Load the path
clear
pth=uigetfile_n_dir();  %OE40 %Final Classfied folder
for i=1:size(pth,2)
    sp=split(pth{1,i},'\');
    [fnmt,ptht]=uigetfile('.mat',['Select the tracking file' sp{end,1}],'Multiselect','on');
fnmtr{i,1}=fnmt;
fnmtr{i,2}=ptht;
end
load('C:\Users\Administrator\Dropbox\In vivo imaging\Analysis\Blind_Answer_total.mat')
save_folder='H:\Image_data\OlympusTPM\Invivo_images_CFC\';
d_rad=60; % Distance of the lipofuscin from cell centroid (the limit to consider to evaluate).
rng=25; % Pixel range to evaluate the mean value of Green signal.
%% Extracting the lipofuscin tracks from Good and Bad mouse. 
clear lip_tr_cre lip_tr_cre_pow
ext_rng=0;
for i=1:size(pth,2)
for k=1:length(fnmtr{i,1})  % Extract credible tracks from good cells condition.
       lip_tr=importdata([fnmtr{i,2} fnmtr{i,1}{1,k}]);
       g=1;
       for ptl=1:max(lip_tr.trackingData(:,1))
           if size(find(lip_tr.trackingData(:,1)==ptl),1)>3 % Extract the track longer than 3.
% && max(lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),6)...
%               .*lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),8).*lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),9))>2000
           tr=lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),:);    
           [Amp ,Amp_arg]=max(tr(:,6).*tr(:,9).*tr(:,8));  % Among the track use the position where the biggest intensity (Amp*sigmaX*sigmaY)
           if Amp_arg+ext_rng>size(tr,1) || Amp_arg-ext_rng<1
           lip_tr_cre{i,k}(g,:)=tr(Amp_arg,:);
           g=g+1;
           else
           lip_tr_cre{i,k}(g:g+ext_rng*2,:)=tr(Amp_arg-ext_rng:Amp_arg+ext_rng,:);
           g=g+ext_rng*2+1;
           end
           end
       end
end
end

%% Load the individual cell images corresponding to the classification.

load(['C:\Users\Administrator\Dropbox\In vivo imaging\Analysis\Blind_Answer.mat'])
for p=1:size(pth,2)
   clear Lip_list MAX_TXN_im Postxn TXN_im
   stem=split(pth{1,p},'\');
   mousenumb{1,p}=char(stem{end,1});
   Cell=importdata(['C:\Users\Administrator\Dropbox\In vivo imaging\Analysis\New\Total\' mousenumb{1,p} '.mat']);
   for j=1:size(BlindAnswer,1)
        if ~isempty(strfind(mousenumb{1,p},char(BlindAnswer(j,1)))) && ~isempty(strfind(char(BlindAnswer(j,1)),mousenumb{1,p}))
            ref=j; %row number that matches to the mouse number
        end
   end
   for j=1:size(fnmtr{p,1},2) %HC CFC Ret days of the image.
   Postxn=[1:1:size(Cell.Data,1)]'; %Regard the cells with TXN at least once.
   clear ds
   for ii=1:size(lip_tr_cre{p,j},1) % Evaluate the distance of lipofuscin near the cell.
       for jj=1:size(Postxn,1) %Filtered cell 
   ds(ii,jj)=distance_BH(lip_tr_cre{p,j}(ii,[3 4 2]),Cell.Data(Postxn(jj,1),1:3));
   end
   end
   imgs=imageDatastore([pth{1,p} '\' char(BlindAnswer(ref,j+1))],'IncludeSubfolders',true,'LabelSource','foldernames'); % Load the list of cropped images
   sp=split(imgs.Files,'\'); sp=split(sp(:,end),'.'); sp=cellfun(@str2num,sp(:,1)); 
   for c=1:size(Postxn,1)  % Cell % Load the z-stack and lipofuscin list.
   iminfo=imfinfo(imgs.Files{find(sp==Postxn(c,1)),1});
       for k=1:numel(iminfo) %TXN zstack
       TXN_im(:,:,k)=imread(imgs.Files{find(sp==Postxn(c,1)),1},k);
       end
   MAX_TXN_im{j}(:,:,c)=max(TXN_im,[],3);
   Lip_list{j,c}=lip_tr_cre{p,j}(find(ds(:,c)<d_rad),:);
   end
   end 

Class_dat=[];
Class_cat=categorical([]);
maxint=[];
clear Lip_int
g=1;
    for j=1:size(MAX_TXN_im,2) %Day
siz=size(MAX_TXN_im{j});
maxint(:,:,g:size(MAX_TXN_im{j},3)+g-1)=MAX_TXN_im{j}(round(siz(1,1)/2)-rng:round(siz(1,1)/2)+rng,round(siz(1,2)/2)-rng:round(siz(1,2)/2)+rng,:); %crop
tmp=cellfun(@(x) x(:,6).*x(:,8).*x(:,9) ,Lip_list(j,:),'un',0);
n_lip=reshape(cell2mat(cellfun(@size, tmp(1,:),'un',0)),2,size(Postxn,1));
n_lip_norm=reshape(cell2mat(cellfun(@size, tmp(1,:),'un',0)),2,size(Postxn,1))./mean(squeeze(mean(mean(maxint(:,:,g:size(MAX_TXN_im{j},3)+g-1),1),2)));
mean_std_lip=[cell2mat(cellfun(@mean, tmp(1,:),'un',0))' cell2mat(cellfun(@sum, tmp(1,:),'un',0))' cell2mat(cellfun(@std, tmp(1,:),'un',0))'];
mean_std_lip_norm=mean_std_lip./mean(squeeze(mean(mean(maxint(:,:,g:size(MAX_TXN_im{j},3)+g-1),1),2)));
% tmp=cellfun(@(x) x(:,6) ,Lip_list(j,:),'un',0);
out=cellfun(@sort,tmp(1,:),'un',0);
 clear top_ten top_ten_norm
 top_ten=zeros(size(tmp,2),1);
for i=1:size(tmp,2)
    if ~isempty(out{1,i}) 
top_ten(i,1) = mean(out{1,i}(round(size(out{1,i},1)*0.9):end));  
    end
end
top_ten_norm=top_ten./mean(squeeze(mean(mean(maxint(:,:,g:size(MAX_TXN_im{j},3)+g-1),1),2)));
Lip_int(g:size(MAX_TXN_im{j},3)+g-1,:)=[n_lip(1,:)' n_lip_norm(1,:)' mean_std_lip mean_std_lip(:,1)./mean_std_lip(:,3) mean_std_lip_norm top_ten top_ten_norm];
g=g+size(MAX_TXN_im{j},3);
    end
    
maxint=double(reshape(maxint,size(maxint,1)*size(maxint,2),1,size(maxint,3)));
Class_dat=[squeeze([mean(maxint) std(maxint) mean(maxint)./std(maxint) ])' Lip_int];
%Class_dat=(Class_dat,size(Cell.Data,1),size(MAX_TXN_im,2));
save([save_folder mousenumb{1,p}],'Class_dat','Lip_list')

thres=[24 4 9.7];
f=[4 10 14];
%plot_3d_thr_plane(Class_dat,f,thres)
clear thres_log
thres_log=zeros(size(Class_dat,1),1);
for i=1:size(f,2)
thres_log=thres_log+double(Class_dat(:,f(1,i))>thres(1,i));
end
thres_log=thres_log==size(f,2);
Cell.thres=reshape(thres_log,size(Cell.Data,1),size(MAX_TXN_im,2));
save([save_folder mousenumb{1,p}],'Class_dat','Lip_list','Cell')
p
end
%%
for p=1%:size(pth,2)
    stem=split(pth{1,p},'\');
   mousenumb{1,p}=char(stem{end,1});
load([save_folder mousenumb{1,p} '.mat'])
thres=[7 4.5 14]; 
thres_HC=[5 3.7 11];
f=[4 10 14]; %Lipofuscinc number, Lip_int/G_mean, top 10% Lip_int/G_mean
behave=2;
clear thres_log
thres_log=zeros(size(Class_dat,1),1);
for i=1:size(f,2)
thres_log(1:size(Cell.Data,1))=thres_log(1:size(Cell.Data,1))+double(Class_dat(1:size(Cell.Data,1),f(1,i))>thres_HC(1,i));
thres_log(size(Cell.Data,1)+1:end)=thres_log(size(Cell.Data,1)+1:end)+double(Class_dat(size(Cell.Data,1)+1:end,f(1,i))>thres(1,i));
end
thres_log=thres_log==size(f,2);
Cell.thres=reshape(thres_log,size(Cell.Data,1),size(fnmtr{p,1},2));

%subplot(3,12,3*(p-1)+1)
subplot(1,3,1)
for i=1:size(fnmtr{p,1},2)
    if i==1
plot_3d_thr_plane(Class_dat((i-1)*size(Cell.Data,1)+1:(i)*size(Cell.Data,1),:),f,thres_HC,0)
    else
hold all
plot_3d_thr_plane(Class_dat((i-1)*size(Cell.Data,1)+1:(i)*size(Cell.Data,1),:),f,thres,0)
    end
end
clear predict imm
%subplot(3,12,3*(p-1)+2)
subplot(1,3,2)
im_fnm=['H:\Image_data\OlympusTPM\Invivo_images_CFC\' mousenumb{1,p} '\Transformed2\crop\' mousenumb{1,p} '_CFC.tif'];
iminfo=imfinfo(im_fnm);
for i=1:numel(iminfo)
imm(:,:,i)=imread(im_fnm,i);
end
imagesc(max(imm,[],3))
axis equal tight off
colormap('gray')
hold all
% plot(PosCell.Data(Postxn(find(sum(predict{2}==1,2)<1),1),1),PosCell.Data(Postxn(find(sum(predict{2}==1,2)<1),1),2),'ro')
% plot(PosCell.Data(Postxn(find(sum(predict{2}==1,2)>1),1),1),PosCell.Data(Postxn(find(sum(predict{2}==1,2)>1),1),2),'bo')
% plot(NegCell.Data(Negtxn(find(sum(predict{2}==1,2)<1),1),1),NegCell.Data(Negtxn(find(sum(predict{2}==1,2)<1),1),2),'ro')
% plot(NegCell.Data(Negtxn(find(sum(predict{2}==1,2)>1),1),1),NegCell.Data(Negtxn(find(sum(predict{2}==1,2)>1),1),2),'bo')
plot(Cell.Data(find(Cell.thres(:,behave)==1),1),Cell.Data(find(Cell.thres(:,behave)==1),2),'ro')
plot(Cell.Data(find(Cell.thres(:,behave)==0),1),Cell.Data(find(Cell.thres(:,behave)==0),2),'bo')
%save('E:\BACKUP\대학원\연구실\MY_Projects\In_vivo_imaging\Cell_detection\TXN_detection_deep_learning\Cell_acceptance\Class_data','Class_table')
%subplot(3,12,3*(p-1)+3)
subplot(1,3,3)
imagesc(Cell.thres)
title(mousenumb{1,p})
name=mousenumb{1,p};
save([save_folder mousenumb{1,p}],'Class_dat','Lip_list','Cell','thres','thres_HC','name')
end