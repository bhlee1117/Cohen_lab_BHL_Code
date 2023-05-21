clear
[fnm,pth]=uigetfile('*.tif','Select the raw image','Multiselect','on');
channel=2;
seg_y=15;
seg_stack=3;
output_filename='spine';
output_channel= [1 4]; % 830 RYC, 960 RYC, 1030 RYC
%%
if ischar(fnm)
    fnm={fnm};
end
ref=1;
for i=1:length(fnm)
imginf=imfinfo([pth char(fnm{1,i})]);
numstack=size(imginf,1);  
sp_name=split(fnm{1,i},'_');
for j=1:numstack
im{i}(:,:,j)=imread([pth char(fnm{1,i})],j);
end
end
%%
figure(1)
for i=1:length(fnm)
subplot(1,length(fnm),i)
imagesc(im{i}(:,:,round(numstack/channel/2)))
hold all
axis equal
end
for i=1:length(fnm)
    subplot(1,length(fnm),i)
    axis equal
   pos(i,:)=ginput(1);
end
close(figure(1))
%%
n=5;
clear ptl_br

for i=1:length(fnm)
ptl_br(i,:)=mean(mean(im{i}(pos(i,2)-n:pos(i,2)+n,pos(i,1)-n:pos(i,1)+n,1:numstack/channel),1),2);
stack_shifted(1,:)=[1 numstack/channel];
if i~=ref
    [cor cor_bin]=crosscorr(ptl_br(ref,:),ptl_br(i,:),'NumLags',numstack/channel-1);
    [m argm]=max(cor);
    shift(i)=cor_bin(argm);
stack_shifted(i,:)=[1-shift(i) numstack/channel-shift(i)];    
end
subplot(1,2,1)
plot([1:1:numstack/channel],ptl_br(i,:))
hold all
subplot(1,2,2)
plot([stack_shifted(i,1):1:stack_shifted(i,2)],ptl_br(i,:))
hold all
end
stack_range=[max(stack_shifted(:,1)) min(stack_shifted(:,2))];



%%
clear im_seg
for i=1:length(fnm)
    for j=1:ceil(size(im{i},1)/seg_y)
    for k=1:ceil((stack_range(1,2)-stack_range(1,1)+1)/seg_stack)
        for c=1:channel
        if j==ceil(size(im{i},1)/seg_y)    
     if k==ceil((stack_range(1,2)-stack_range(1,1)+1)/seg_stack)
    im_seg{i,c,j,k}=im{i}(seg_y*(j-1)+1:size(im{i},1),:,numstack/channel*(c-1)+seg_stack*(k-1)+1+shift(i):numstack/channel*(c-1)+stack_range(1,2)+shift(i));                 
     else
    im_seg{i,c,j,k}=im{i}(seg_y*(j-1)+1:size(im{i},1),:,numstack/channel*(c-1)+seg_stack*(k-1)+stack_range(1,1)+shift(i):numstack/channel*(c-1)+stack_range(1,1)-1+seg_stack*k+shift(i));        
     end
     else
         if  k==ceil((stack_range(1,2)-stack_range(1,1)+1)/seg_stack)
    im_seg{i,c,j,k}=im{i}(seg_y*(j-1)+1:seg_y*j,:,numstack/channel*(c-1)+seg_stack*(k-1)+1+shift(i):numstack/channel*(c-1)+stack_range(1,2)+shift(i));
         else
             
    im_seg{i,c,j,k}=im{i}(seg_y*(j-1)+1:seg_y*j,:,numstack/channel*(c-1)+seg_stack*(k-1)+stack_range(1,1)+shift(i):numstack/channel*(c-1)+stack_range(1,1)-1+seg_stack*k+shift(i));
     end
        end
        end
    end
    end
end
%%
clear im_rec
h=waitbar(0,'Reconstruction');
for i=1:size(im_seg,1)
    
    im_rec{i}=zeros(size(im{i},1),size(im{i},2),size(im_seg,4));
    for j=1:size(im_seg,3)
    for k=1:size(im_seg,4)
         clear tdcorr
    for z=1:size(im_seg{i,1,j,k},3) %{:,c,:,:} -> set c as reference channel.
        ave_im=mean(im_seg{i,1,j,k},3);
       %tdcorr(z)=max(max(xcorr2(ave_im,im_seg{i,j,k}(:,:,z))));
       xcf=xcorr2(ave_im,im_seg{i,1,j,k}(:,:,z));
        if j~=size(im_seg,3)
       tdcorr(z)=xcf(seg_y,size(im_seg{i,1,j,k},2));
        else
       tdcorr(z)=xcf(size(ave_im,1),size(im_seg{i,1,j,k},2));
        end
    end
    [m maxz]=max(tdcorr);
    for c=1:channel
        if j~=size(im_seg,3)
    im_rec{i}(seg_y*(j-1)+1:seg_y*j,:,k+size(im_seg,4)*(c-1))=im_seg{i,c,j,k}(:,:,maxz);
        else
    im_rec{i}(1:size(ave_im,1),:,k+size(im_seg,4)*(c-1))=im_seg{i,c,j,k}(:,:,maxz);
        end
    end
    end
    waitbar((1/size(im_seg,1)*(i-1)+1/size(im_seg,1)*j/size(im_seg,3)))
    end
    
    for z=1:size(im_rec{i},3)
imwrite(uint16(im_rec{i}(:,:,z)),[pth,'Corr_',fnm{1,i}],'Writemode','Append','Compression','none');
    end
end
close(h)

%%
clear im_rec_trs
% for i=1:size(im_rec{ref},3)
%     imwrite(uint16(im_rec{ref}(:,:,i)),[pth,output_filename,'.tif'],'Writemode','Append','Compression','none');
% end
for i=1:size(im_rec{ref},3)/channel
    for c=1:channel
    im_rec_trs{c}(:,:,i)=im_rec{ref}(:,:,i+size(im_rec{ref},3)/channel*(c-1));
    end
end
trans_stack=5;

i=1;
while i<=size(fnm,2)
    if i~=ref
% BaseFile = im_rec{ref}(:,:,trans_stack);
% FileToWarp = im_rec{i}(:,:,trans_stack);
% % FileToWarp_lip = [PathName1 fnm{1,i}(5:end)];

%WarpImageBaseFile = sscanf(FileToWarp,'%c',size(FileToWarp,2)-4);
%imginf=imfinfo(FileToWarp);

ImageToWarp = im_rec{i}(:,:,trans_stack);
BaseImage = im_rec{ref}(:,:,trans_stack);

tempWarp = double(ImageToWarp);
tempBase = double(BaseImage);

tempWarp = tempWarp./max(max(tempWarp))*3;
tempBase = tempBase./max(max(tempBase))*3;

% Select the control points.
h = cpselect(tempWarp, tempBase);
uiwait(msgbox('Click OK after closing the CPSELECT window.','Waiting...'))

movingPoints = cpcorr(movingPoints, fixedPoints,tempWarp,tempBase);
tform2 = maketform('projective',movingPoints, fixedPoints);

for z=1:size(im_rec{i},3)
    if mod(z-1,size(im_rec{i},3)/channel)==0
        cc=size(im_rec_trs,2)+1;
    end

    ImageToWarp = im_rec{i}(:,:,z);
[transVis xdata ydata] = imtransform(ImageToWarp, tform2);
transVisWithOffset = imtransform(ImageToWarp, tform2, 'XData', [1 (size(ImageToWarp,2)+tform2.tdata.T(3,1))],'YData', [1 (size(ImageToWarp,1)+ tform2.tdata.T(3,2))]);
transVisBuffrd = uint16(zeros(size(BaseImage)));
transVisBuffrd(1:size(transVisWithOffset,1),1:size(transVisWithOffset,2)) = transVisWithOffset;
Transformation_Matrix=tform2.tdata.T(1:2,1:2);
Translation_Vector=tform2.tdata.T(3,1:2)';
transVisBuffrd = transVisBuffrd((1:size(ImageToWarp,1)), (1:size(ImageToWarp,2)));

c=ceil(z/(size(im_rec{i},3)/channel));

im_rec_trs{cc}(:,:,z-(c-1)*(size(im_rec{i},3)/channel))=transVisBuffrd;
%imwrite(uint16(transVisBuffrd),[pth,output_filename,'.tif'],'Writemode','Append','Compression','none');
%imwrite(transVisBuffrd,[save_path 'T_' char(sp_name{6,1})],'tif','Writemode','Append','Compression', 'none');

% ImageToWarp = imread(FileToWarp_lip,'tif',z);
% [transVis xdata ydata] = imtransform(ImageToWarp, tform2);
% transVisWithOffset = imtransform(ImageToWarp, tform2, 'XData', [1 (size(ImageToWarp,2)+tform2.tdata.T(3,1))],'YData', [1 (size(ImageToWarp,1)+ tform2.tdata.T(3,2))]);
% transVisBuffrd = uint16(zeros(size(BaseImage)));
% transVisBuffrd(1:size(transVisWithOffset,1),1:size(transVisWithOffset,2)) = transVisWithOffset;
% Transformation_Matrix=tform2.tdata.T(1:2,1:2);
% Translation_Vector=tform2.tdata.T(3,1:2)';
% transVisBuffrd = transVisBuffrd((1:size(ImageToWarp,1)), (1:size(ImageToWarp,2)));
% 
% imwrite(transVisBuffrd,[save_path 'T_lip_' char(sp_name{6,1})],'tif','Writemode','Append','Compression', 'none');
 
end

%save([WarpImageBaseFile,'T.mat'], 'tform2');
    end
    i=i+1;
end
%% 
for i=1:size(im_rec_trs,2)
    sz(i,:)=[size(im_rec_trs{i},1) size(im_rec_trs{i},2)];
end
min_sz=min(sz);
for i=output_channel
    
for z=1:size(im_rec_trs{i},3)
imwrite(uint16(imcrop(im_rec_trs{i}(:,:,z),[1 1 min_sz(1,1) min_sz(1,2)])),[pth,output_filename,'.tif'],'tif','Writemode','Append','Compression', 'none');
end
end
save([pth,output_filename,'.mat'],'im_rec_trs')