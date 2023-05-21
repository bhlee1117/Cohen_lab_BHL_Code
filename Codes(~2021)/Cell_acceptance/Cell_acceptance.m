%% Load the path
clear
pth=uigetfile_n_dir();  %OE40 %Final Classfied folder
[fnmtr,pthtr]=uigetfile('.mat','Select the tracking file','Multiselect','on');
pth_neg=uigetfile_n_dir(); %OE44
[fnmntr,pthntr]=uigetfile('.mat','Select the negative tracking file','Multiselect','on');
PosCell=importdata('C:\Users\Administrator\Dropbox\In vivo imaging\Analysis\New\Total\OE39.mat');
NegCell=importdata('C:\Users\Administrator\Dropbox\In vivo imaging\Analysis\New\Total\ROE27.mat');
load('C:\Users\Administrator\Dropbox\In vivo imaging\Analysis\Blind_Answer_total.mat')
d_rad=60; % Distance of the lipofuscin from cell centroid (the limit to consider to evaluate).
%% Extracting the lipofuscin tracks from Good and Bad mouse. 
clear lip_tr_cre lip_tr_cre_pow
ext_rng=0;
for k=1:length(fnmtr)  % Extract credible tracks from good cells condition.
       lip_tr=importdata([pthtr fnmtr{1,k}]);
       g=1;
       for ptl=1:max(lip_tr.trackingData(:,1))
           if size(find(lip_tr.trackingData(:,1)==ptl),1)>3 % Extract the track longer than 3.
% && max(lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),6)...
%               .*lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),8).*lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),9))>2000
           tr=lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),:);    
           [Amp ,Amp_arg]=max(tr(:,6).*tr(:,9).*tr(:,8));  % Among the track use the position where the biggest intensity (Amp*sigmaX*sigmaY)
           if Amp_arg+ext_rng>size(tr,1) || Amp_arg-ext_rng<1
           lip_tr_cre{1,k}(g,:)=tr(Amp_arg,:);
           g=g+1;
           else
           lip_tr_cre{1,k}(g:g+ext_rng*2,:)=tr(Amp_arg-ext_rng:Amp_arg+ext_rng,:);
           g=g+ext_rng*2+1;
           end
           end
       end
end
for k=1:length(fnmntr) % Extract credible tracks from bad cells condition.
       lip_tr=importdata([pthntr fnmntr{1,k}]);
       g=1;
       for ptl=1:max(lip_tr.trackingData(:,1))
           if size(find(lip_tr.trackingData(:,1)==ptl),1)>3 % Extract the track longer than 3.
%                && max(lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),6)...
%              .*lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),8).*lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),9))>2000
           tr=lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),:);    
           [Amp ,Amp_arg]=max(tr(:,6).*tr(:,9).*tr(:,8));
           if Amp_arg+ext_rng>size(tr,1) || Amp_arg-ext_rng<1
           lip_tr_cre{2,k}(g,:)=tr(Amp_arg,:);
           g=g+1;
           else
           lip_tr_cre{2,k}(g:g+ext_rng*2,:)=tr(Amp_arg-ext_rng:Amp_arg+ext_rng,:);
           g=g+ext_rng*2+1;
           end
       end
       end
end
%% Load the individual cell images corresponding to the classification.

clear MAX_TXN_im MAX_TXN_im Lip_list_TXN Lip_list_NXN MAX_NXN_im Postxn
load(['C:\Users\Administrator\Dropbox\In vivo imaging\Analysis\Blind_Answer.mat'])
    %BlindAnswer=table2array(BlindAnswer);

   stem=split(pth{1,1},'\');
   mousenumb{1,1}=char(stem{end,1});
   for j=1:size(BlindAnswer,1)
        if ~isempty(strfind(mousenumb{1,1},char(BlindAnswer(j,1)))) && ~isempty(strfind(char(BlindAnswer(j,1)),mousenumb{1,1}))
            ref=j; %row number that matches to the mouse number
        end
   end
   for j=1:size(fnmtr,2) %HC CFC Ret days of the image.
   Postxn=find(sum(PosCell.Data(:,4:end)==2,2)>0); %Regard the cells with TXN at least once.
   clear ds
   for ii=1:size(lip_tr_cre{1,j},1) % Evaluate the distance of lipofuscin near the cell.
       for jj=1:size(Postxn,1) %Filtered cell 
   ds(ii,jj)=distance_BH(lip_tr_cre{1,j}(ii,[3 4 2]),PosCell.Data(Postxn(jj,1),1:3));
   end
   end
   
   imgs=imageDatastore([pth{1,1} '\' char(BlindAnswer(ref,j+1))],'IncludeSubfolders',true,'LabelSource','foldernames'); % Load the list of cropped images
   sp=split(imgs.Files,'\'); sp=split(sp(:,end),'.'); sp=cellfun(@str2num,sp(:,1)); 
   for c=1:size(Postxn,1)  % Cell % Load the z-stack and lipofuscin list.
   iminfo=imfinfo(imgs.Files{find(sp==Postxn(c,1)),1});
       for k=1:numel(iminfo) %TXN zstack
       TXN_im(:,:,k)=imread(imgs.Files{find(sp==Postxn(c,1)),1},k);
       end
   MAX_TXN_im{j}(:,:,c)=max(TXN_im,[],3);
   Lip_list_TXN{j,c}=lip_tr_cre{1,j}(find(ds(:,c)<d_rad),:);
   end
   
   end 
 %%%%%% Same for negative population
 
  stem=split(pth_neg{1,1},'\');
   mousenumb{1,1}=char(stem{end,1});
    for j=1:size(BlindAnswer,1)
        if ~isempty(strfind(mousenumb{1,1},char(BlindAnswer(j,1)))) && ~isempty(strfind(char(BlindAnswer(j,1)),mousenumb{1,1}))
            ref=j; %row number that matches to the mouse number
        end
    end
   for j=1:size(fnmntr,2) %HC CFC Ret days of the image.
   %Negtxn=find(sum(NegCell.Data(:,4:end)==1,2)==size(NegCell.Data,2)-3); %Load the cell that never has TXN site.
   Negtxn=[1:1:size(NegCell.Data,1)]'; % All cells 
   clear ds
   for ii=1:size(lip_tr_cre{2,j},1) % Evaluate the distance of lipofuscin near the cell.
       for jj=1:size(Negtxn,1)
   ds(ii,jj)=distance_BH(lip_tr_cre{2,j}(ii,[3 4 2]),NegCell.Data(Negtxn(jj,1),1:3));
   end
   end
   imgs=imageDatastore([pth_neg{1,1} '\' char(BlindAnswer(ref,j+1))],'IncludeSubfolders',true,'LabelSource','foldernames'); % Load the cells rearranging the blind.
   sp=split(imgs.Files,'\'); sp=split(sp(:,end),'.'); sp=cellfun(@str2num,sp(:,1));
   for c=1:size(Negtxn,1) %Load the bad cell and lipofuscin list.
   iminfo=imfinfo(imgs.Files{find(sp==Negtxn(c,1)),1});
       for k=1:numel(iminfo) %TXN zstack
       NXN_im(:,:,k)=imread(imgs.Files{find(sp==Negtxn(c,1)),1},k);
       end
   MAX_NXN_im{j}(:,:,c)=max(NXN_im,[],3);
   Lip_list_NXN{j,c}=lip_tr_cre{2,j}(find(ds(:,c)<d_rad),:);
   end
   end
%% Extract and make a training data [Lipo list, Green intensity]
rng=25;
Class_dat=[];
Class_cat=categorical([]);
maxint=[];
clear Lip_int
%Lip_int=[];
g=1;
    for j=1:size(MAX_TXN_im,2) %Day
siz=size(MAX_TXN_im{j});
maxint(:,:,g:size(MAX_TXN_im{j},3)+g-1)=MAX_TXN_im{j}(round(siz(1,1)/2)-rng:round(siz(1,1)/2)+rng,round(siz(1,2)/2)-rng:round(siz(1,2)/2)+rng,:); %crop
tmp=cellfun(@(x) x(:,6).*x(:,8).*x(:,9) ,Lip_list_TXN(j,:),'un',0);
n_lip=reshape(cell2mat(cellfun(@size, tmp(1,:),'un',0)),2,size(Postxn,1));
n_lip_norm=reshape(cell2mat(cellfun(@size, tmp(1,:),'un',0)),2,size(Postxn,1))./mean(squeeze(mean(mean(maxint(:,:,g:size(MAX_TXN_im{j},3)+g-1),1),2)));
mean_std_lip=[cell2mat(cellfun(@mean, tmp(1,:),'un',0))' cell2mat(cellfun(@sum, tmp(1,:),'un',0))' cell2mat(cellfun(@std, tmp(1,:),'un',0))'];
mean_std_lip_norm=mean_std_lip./mean(squeeze(mean(mean(maxint(:,:,g:size(MAX_TXN_im{j},3)+g-1),1),2)));
% tmp=cellfun(@(x) x(:,6) ,Lip_list_TXN(j,:),'un',0);
out=cellfun(@sort,tmp(1,:),'un',0);
 clear top_ten
for i=1:size(tmp,2) 
top_ten(i,1) = mean(out{1,i}(round(size(out{1,i},1)*0.9):end));  
end
top_ten_norm=top_ten./mean(squeeze(mean(mean(maxint(:,:,g:size(MAX_TXN_im{j},3)+g-1),1),2)));
clear ttmp ttmp_norm
% for i=1:size(tmp,2) 
%     ttmp(i,:)=histcounts(tmp{1,i},[0:2000:10000],'Normalization','Probability'); 
%     ttmp_norm(i,:)=histcounts(tmp{1,i}./mean(mean(maxint(:,:,(j-1)*size(tmp,2)+i))),[0:4:20],'Normalization','Probability');
% end
mean(squeeze(mean(mean(maxint(:,:,g:size(MAX_TXN_im{j},3)+g-1),1),2)))
Lip_int(g:size(MAX_TXN_im{j},3)+g-1,:)=[n_lip(1,:)' n_lip_norm(1,:)' mean_std_lip mean_std_lip(:,1)./mean_std_lip(:,3) mean_std_lip_norm top_ten top_ten_norm];
%Lip_int(g:size(MAX_TXN_im{j},3)+g-1,1)=cell2mat(cellfun(@mean,cellfun(@(x) x(:,6).*x(:,8).*x(:,9) ,Lip_list_TXN(j,:),'un',0),'un',0))'; %intensity mean
Class_cat(g:size(MAX_TXN_im{j},3)+g-1,1)=categorical(zeros(size(MAX_TXN_im{j},3),1)+1,1,{'Pos'});
g=g+size(MAX_TXN_im{j},3);
    end
    for j=1:size(MAX_NXN_im,2) %Day
        siz=size(MAX_NXN_im{j});
maxint(:,:,g:size(MAX_NXN_im{j},3)+g-1)=MAX_NXN_im{j}(round(siz(1,1)/2)-rng:round(siz(1,1)/2)+rng,round(siz(1,2)/2)-rng:round(siz(1,2)/2)+rng,:); 
tmp=cellfun(@(x) x(:,6).*x(:,8).*x(:,9) ,Lip_list_NXN(j,:),'un',0);
n_lip=reshape(cell2mat(cellfun(@size, tmp(1,:),'un',0)),2,size(Negtxn,1));
n_lip_norm=reshape(cell2mat(cellfun(@size, tmp(1,:),'un',0)),2,size(Negtxn,1))./mean(squeeze(mean(mean(maxint(:,:,g:size(MAX_NXN_im{j},3)+g-1),1),2)));
mean_std_lip=[cell2mat(cellfun(@mean, tmp(1,:),'un',0))' cell2mat(cellfun(@sum, tmp(1,:),'un',0))' cell2mat(cellfun(@std, tmp(1,:),'un',0))'];
mean_std_lip_norm=mean_std_lip./mean(squeeze(mean(mean(maxint(:,:,g:size(MAX_NXN_im{j},3)+g-1),1),2)));
out=cellfun(@sort,tmp(1,:),'un',0);
 clear top_ten
for i=1:size(tmp,2) 
    if ~isempty(out{1,i})
top_ten(i,1) = mean(out{1,i}(round(size(out{1,i},1)*0.9):end));  
    end
end
top_ten_norm=top_ten./mean(squeeze(mean(mean(maxint(:,:,g:size(MAX_NXN_im{j},3)+g-1),1),2)));
%tmp=cellfun(@(x) x(:,6) ,Lip_list_NXN(j,:),'un',0);
% 
clear ttmp ttmp_norm
% for i=1:size(tmp,2) 
%     ttmp(i,:)=histcounts(tmp{1,i},[0:2000:10000],'Normalization','Probability'); 
%     ttmp_norm(i,:)=histcounts(tmp{1,i}./mean(mean(maxint(:,:,(j-1)*size(tmp,2)+i))),[0:4:20],'Normalization','Probability');
% end
Lip_int(g:size(MAX_NXN_im{j},3)+g-1,:)=[n_lip(1,:)' n_lip_norm(1,:)' mean_std_lip mean_std_lip(:,1)./mean_std_lip(:,3) mean_std_lip_norm top_ten top_ten_norm];
Class_cat(g:size(MAX_NXN_im{j},3)+g-1,1)=categorical(zeros(size(MAX_NXN_im{j},3),1)+1,1,{'Neg'});
g=g+size(MAX_NXN_im{j},3);
end
    
maxint=double(reshape(maxint,size(maxint,1)*size(maxint,2),1,size(maxint,3)));

Class_dat=[squeeze([mean(maxint) std(maxint) mean(maxint)./std(maxint) ])' Lip_int];
%Class_dat(size(Class_dat,1)+1:size(Class_dat,1)+size(maxint,3),:)=squeeze(hist(maxint,[50:50:2500]))';

Class_table=table(Class_dat,Class_cat);

save('E:\BACKUP\대학원\연구실\MY_Projects\In_vivo_imaging\Cell_detection\TXN_detection_deep_learning\Cell_acceptance\Class_data','Class_table')
%%
thres=[20 4 9.7];
f=[4 10 14];
plot_3d_thr_plane(Class_dat,f,thres,0)

%%
thres_log=zeros(size(Class_dat,1),1);
for i=1:size(f,2)
thres_log=thres_log+double(Class_dat(:,f(1,i))>thres(1,i));
end
thres_log=thres_log==size(f,2);
%%
% [trainedClassifier, validationAccuracy] = trainClassifier(Class_table);
% yfit = trainedClassifier.predictFcn(Class_table);
clear predict imm
for i=1:63
imm(:,:,i)=imread('H:\Image_data\OlympusTPM\Invivo_images_CFC\ROE27\Transformed2\crop\ROE27_CFC.tif',i);
end
imagesc(max(imm,[],3))
axis equal tight off
colormap('gray')
hold all
predict{1}=reshape(thres_log(1:size(Lip_list_TXN,1)*size(Lip_list_TXN,2),1),size(Lip_list_TXN,2),size(Lip_list_TXN,1));
predict{2}=reshape(thres_log(size(Lip_list_TXN,1)*size(Lip_list_TXN,2)+1:end,1),size(Lip_list_NXN,2),size(Lip_list_NXN,1));
% plot(PosCell.Data(Postxn(find(sum(predict{2}==1,2)<1),1),1),PosCell.Data(Postxn(find(sum(predict{2}==1,2)<1),1),2),'ro')
% plot(PosCell.Data(Postxn(find(sum(predict{2}==1,2)>1),1),1),PosCell.Data(Postxn(find(sum(predict{2}==1,2)>1),1),2),'bo')
% plot(NegCell.Data(Negtxn(find(sum(predict{2}==1,2)<1),1),1),NegCell.Data(Negtxn(find(sum(predict{2}==1,2)<1),1),2),'ro')
% plot(NegCell.Data(Negtxn(find(sum(predict{2}==1,2)>1),1),1),NegCell.Data(Negtxn(find(sum(predict{2}==1,2)>1),1),2),'bo')
plot(NegCell.Data(Negtxn(find(predict{2}(:,2)==1),1),1),NegCell.Data(Negtxn(find(predict{2}(:,2)==1),1),2),'ro')
plot(NegCell.Data(Negtxn(find(predict{2}(:,2)==0),1),1),NegCell.Data(Negtxn(find(predict{2}(:,2)==0),1),2),'bo')
% plot(PosCell.Data(Postxn(find(predict{1}(:,4)==1),1),1),PosCell.Data(Postxn(find(predict{1}(:,4)==1),1),2),'ro')
% plot(PosCell.Data(Postxn(find(predict{1}(:,4)==0),1),1),PosCell.Data(Postxn(find(predict{1}(:,4)==0),1),2),'bo')
save('E:\BACKUP\대학원\연구실\MY_Projects\In_vivo_imaging\Cell_detection\TXN_detection_deep_learning\Cell_acceptance\Class_data','Class_table')
