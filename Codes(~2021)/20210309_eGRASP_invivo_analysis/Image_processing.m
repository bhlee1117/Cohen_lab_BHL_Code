%% Parameters
clear;
N_day=2; % days
ch=4; % channels
pix_xyz=[0.15 1]; % pixel size in um
save_pth=['D:\Byunghun_Lee\eGRASP\Images\'];
%% Load images
for d=1:N_day
    [filename{d}] = startupGUI_GRASP([num2ordinal(d) ' day']);
end
[n t]=cellfun(@size, filename);
if length(unique(n))>1
    error('Check the loaded files!')
else
    N=unique(n);
    for d=1:N_day % day
        for fn=1:N % Areas
            numim=size(imfinfo(filename{d}{fn}),1);
            Zlength=numim/ch;
            ind=reshape([1:numim],ch,[]);
            for i=1:numim
                [a b]=find(ind==i);
                Im{d,fn}(:,:,b,a)=imread(filename{d}{fn},i);
            end
            [ROI cc]=imcrop(mean(Im{d,fn}(:,:,:,3),3)/max(max(mean(Im{d,fn}(:,:,:,3),3))));
            Im{d,fn}=Im{d,fn}(round(cc(2)):round(cc(2))+round(cc(4))-1,round(cc(1)):round(cc(1))+round(cc(3))-1,:,:);
        end
    end
end
close all
% Image alignment btw days
%% Manually select dendrite end coordinates
for d=1%:size(Im,1)
    for fn=1%1:size(Im,2)
        [c_image_xy]=xyz_slice(Im{d,fn});
        D_coor{d,fn}={[]};
        D_coor{d,fn}=Dendrite_selection(Im{d,fn},c_image_xy,pix_xyz,'add',D_coor{d,fn});
    end
end
%save(['D:\Byunghun_Lee\eGRASP\Images\20210310_data.mat'],'D_coor','Im');
%% Find the same dendrites in D+2 and D+6 images
    for fn=1:size(Im,2)
       [candidatePos{fn} Y{fn} reg_angle(fn) shift(fn)]=find_rng_different_image(D_coor{fn},Im(:,fn));
    end
    
    load(['D:\Byunghun_Lee\eGRASP\Images\20210310_data.mat'],'D_coor','Im','reg_angle','shift','Y','candidatePos');
%% Extract dendrite intensity of mBeRFP channel on day -1, +2 and +6
for fn=1:size(Im,2)
Dend_int{fn}=extract_dend_intensity(Y{fn},D_coor{fn},candidatePos{fn},pix_xyz);
M=cellfun(@mean,Dend_int{fn}); S=cellfun(@std,Dend_int{fn}); [N N2]=cellfun(@size,Dend_int{fn});
figure
for i=1:size(Dend_int{fn},2)
errorbar(M(1,i),M(2,i),S(2,i),S(2,i),S(1,i),S(1,i),'marker','o')
hold all
end
end
 save(['D:\Byunghun_Lee\eGRASP\Images\20210310_data.mat'],'D_coor','Im','reg_angle','shift','Y','candidatePos','Y3','Dend_int');
 
 %% show projected image
Y3=projection_z(Y{fn},D_coor{fn},candidatePos{fn},pix_xyz); 
 