% VR-Arc & Calcium imaging, registration and data extraction code
%
% Website : http://neurobiophysics.snu.ac.kr/

% INPUTS 

% Image stacks : 900 nm Arc images 
% NoRMCorr registrated timelapse calcium image

% OUTPUTS
% Cell positions and spatial components 

% MODIFICATION HISTORY : 
%     2020.07.23.
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy, Seoul National University, 2020.07.23.
%% load Arc images (stitched)
clear
[fnmArc pthArc]=uigetfile('*.tif','Select the stitched zstack images','Multiselect','on');
[pthCa]=uigetfile_n_dir();  % Select the .tif file folder taken from the olympus Two photon microscope. 
%[fnmVR pthVR]=uigetfile('*.mat','Select the VR data','Multiselect','on');
save_pth='\\Neurobiophysics\Virtual Reality\VR_image_data\JOE19\';
%%
clear im
for i=1:size(fnmArc,2) % Load Arc images
    im_info=imfinfo([pthArc fnmArc{1,i}]);
    for z=1:numel(im_info)
    im{i}(:,:,z)=imread([pthArc fnmArc{1,i}],z);
    end
end
%% Z stack orientation & registration across multiple days

cmap=distinguishable_colors(size(fnmArc,2));
im_aligned=zstack_align(im,2,30,20,12);  %z_ref,rng,max_shift
merge_image(im_aligned,cmap,fnmArc)
for i=1:size(im_aligned{1},3)
    imwrite(uint16(im_aligned{1}(:,:,i)),[save_pth 'ref_image.tif'],'Writemode','Append')
end
%% Cell detection

cell_list=cell_detection(im_aligned,5,'Add'); %'Add or Find'
%%
cal_reg_done=true;
if ~cal_reg_done
    Calcium_registration(pthCa,2,1,[],0) %pth,ch,ref,step,save_tif 
end
for p=1:size(pthCa,2)
load([pthCa{1,p}(1,1:end-10) '_calcium_image.mat'],'MG','MR')
%MG(MG>60)=0;
tic;[regist.angle(p,1) regist.argz(p,1) regist.roi(p,:)]=find_matZ(MG,im_aligned{1},40);toc;
reg(:,:,1)=im_aligned{1}(:,:,regist.argz(p,1)); 
MG_crop=imcrop(imrotate(imresize(MG,2),regist.angle(p,1),'crop'),regist.roi(p,:)*2);
reg(1:size(MG_crop,1),1:size(MG_crop,2),2)=MG_crop;
merge_image({reg(:,:,1),reg(:,:,2)*3},[1 0 0;1 1 1],{'Arc','Calcium'})
%slider_cell_range(imresize(MG,2),argz(p,1),angle(p,1),2*roi(p,:),cell_list)
end
save([save_pth 'JOE20_Cell.mat'],'cell_list','im_aligned','regist','-v7.3')
%%
load([save_pth 'JOE19_Cell.mat']);
for p=1:size(pthCa,2)
load([pthCa{1,p}(1,1:end-10) '_calcium_image.mat']);
Y=cal_image_rot_crop(Y,regist.angle(p,1),regist.roi(p,:));
result{p}=run_CaImAn(Y,200,1,cell_list,regist.argz(p,1),20); % Calcium image, number of neurons, minSNR, cell_list,argz,rangez
end
for p=1:size(pthCa,2)
    figure
[ref_result{p} baseX{p}]=find_baseline(result{p});
[ref_result{p}.cal_sigma ref_result{p}.cal_transient ref_result{p}.ini_fin]=cal_transient_detection(ref_result{p},2,0.5,baseX{p},1);
end
show_transient(ref_result{1},[50:80],1/30,1)
%%
for i=1:size(pthCa,2)
    load([pthCa{1,i}(1,1:end-10) '_calcium_image.mat'],'MG','MR')
    cal_centroid=cell2mat(cellfun(@mean_BH,ref_result{i}.Coor','UniformOutput',false))';
    MG=cal_image_rot_crop(MG,regist.angle(i,1),regist.roi(i,:));
    [match_CS{i} unmatched_CS{i} match_center{i}]=coordinate_match(cell_list,cal_centroid(:,[1 2]),40,regist.argz(i,1),30,imresize(MG,2)); % max dist, max_Z range
end

[identified_list identified_calcium identified_SpC identified_spike]=identify_sp_comp(ref_result,match_CS);
plot_contours_multipleday(identified_SpC,result{1}.options,[cell_list(:,1:2)/2 cell_list(:,3)],identified_list,regist.argz,[1 2 3])
% plot_place_field(identified_calcium,1,'NormCdf',Position)
save([save_pth 'JOE19_calcium_results.mat'],'result','ref_result','baseX','match_CS','identified_calcium','identified_list','identified_SpC','identified_spike','-v7.3')