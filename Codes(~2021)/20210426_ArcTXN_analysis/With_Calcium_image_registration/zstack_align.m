% z alignment across different days
% Website : http://neurobiophysics.snu.ac.kr/

% INPUTS 

% Image stacks : 900 nm Arc images 
% NoRMCorr registrated timelapse calcium image

% OUTPUTS
% Cell positions and spatial components

% MODIFICATION HISTORY : 
%     2020.07.23.
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy, Seoul National University, 2020.07.23.
function [im_aligned]=zstack_align(im,rr,z_ref,rng,max_shift)
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 700*size(im,2) 700]);
[d1 d2 d3]=cellfun(@size,im);
[m ref]=max(d3);
cmap=[1 0 0;0 1 0;0 0 1;0 1 1];
[maxd]=max([d1' d2']);

for i=fliplr([1:size(im,2)])
I=zeros(maxd(1,1),maxd(1,2),d3(1,i));   
I(1:d1(1,i),1:d2(1,i),1:d3(1,i))=im{i};
im{i}=I;
subplot(1,size(im,2),i)
imagesc(max(im{i}(:,:,z_ref-rng:z_ref+rng),[],3))
colormap('gray');
axis equal tight off
hold all
end

for j=1:4
subplot(1,size(im,2),1)
[crop_im roi_im]=imcrop(max(im{1},[],3)/250); crop_im=crop_im*250;
candidatePos{1}(j,:)=roi_im(:,1:2)+roi_im(:,3:4)/2;
draw_rectangle(roi_im,2,cmap(j,:));
for i=2:size(im,2)
    subplot(1,size(im,2),i)
[ match_data, ~ ] = matchPattern( max(im{i},[],3), crop_im, 0.35,2);  pattRows = size(crop_im,1);  pattCols = size(crop_im,2);
if size(match_data,1)>1
    [m argm]=max(match_data(:,3));
    match_data=match_data(argm,:);
end
match_data(:,1) = match_data(:,1)+ceil(pattRows/2);  match_data(:,2) = match_data(:,2)+ceil(pattCols/2);  candidatePos{i}(j,:) =[match_data(:,2) match_data(:,1)];
draw_rectangle([match_data(:,2)-pattCols/2 match_data(:,1)-pattRows/2 pattCols-1 pattRows-1],2,cmap(j,:));
end
end
pause(0.5)
tic;
f = waitbar(1/size(im,2),'Registration on progress');

for i=2:size(im,2)

movingPoints = cpcorr(candidatePos{i}, candidatePos{1},max(im{i},[],3),max(im{1},[],3));
tform2 = maketform('projective',movingPoints, candidatePos{1});
im{i}=Warp_image(im{i},tform2);
 waitbar(i/size(im,2),f,'Registration on progress');
end
close(f)
toc;
ptl_br{ref}=squeeze(mean(mean(im{ref}(round(candidatePos{1}(1,2)):round(candidatePos{1}(1,2))+10,round(candidatePos{1}(1,1)):round(candidatePos{1}(1,1))+10,:),1),2));
figure2 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 400 400]);
for i=1:size(im,2)
ptl_br{i}   = squeeze(mean(mean(   im{i}(round(candidatePos{rr}(1,2))-20:round(candidatePos{rr}(1,2))+20,round(candidatePos{rr}(1,1))-20:round(candidatePos{rr}(1,1))+20,:),1),2));
stack_shifted(ref,:)=[1 d3(1,i)];
if i~=ref
    [cor cor_bin]=crosscorr(ptl_br{ref},ptl_br{i},'numLags',d3(1,i)-1);
    %[m argm]=max(cor);
    [m argm]=sort(cor,'descend'); 
    g=1;
    while abs(cor_bin(argm(g)))>max_shift
        g=g+1;
    end
    shift(i)=cor_bin(argm(g));
stack_shifted(i,:)=[1-shift(i) d3(1,i)-shift(i)];    
end

subplot(1,2,1)
plot([1:1:d3(1,i)],ptl_br{i})
hold all
subplot(1,2,2)
plot([stack_shifted(i,1):1:stack_shifted(i,2)],ptl_br{i})
hold all
end
[stack_ini ini]=max(stack_shifted(:,1)); [stack_end endd]=min(stack_shifted(:,2));
for i=1:size(im,2)
    im_aligned{i}=double(im{i}(:,:,stack_ini-stack_shifted(i,1)+1:stack_end-stack_shifted(i,1)+1));
end
mat_aligned=cell2mat(im_aligned);
options_rigid = NoRMCorreSetParms('d1',size(im_aligned{1},1),'d2',size(im_aligned{1},2),'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200,'print_msg',0);
f = waitbar(1/size(im,2),'Precise registration on progress');
for i=1:size(im_aligned{1},3)
    waitbar(i/size(im_aligned{1},3),f,['Precise registration on progress ' num2str(i) ' / ' num2str(size(im_aligned{1},3))]);
    [temp_im,shifts1{i},template1,options_rigid] = normcorre(reshape(mat_aligned(:,:,i),size(mat_aligned,1),size(im_aligned{1},2),[]),options_rigid);
    for j=1:size(temp_im,3)
    im_aligned{j}(:,:,i)=temp_im(:,:,j);
    end
end
close(f)

end

