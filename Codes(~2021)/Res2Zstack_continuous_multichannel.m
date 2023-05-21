%%
%     Website : http://neurobiophysics.snu.ac.kr/
% This code makes averaged image from resonant scanning image. 
% Resonant scanned image -> Motion correction (NoRMCorr) -> average -> save in tiff 

% INPUTS 
% Olympus two photon image path, tif format. 
% 

% OUTPUTS
% Averaged image saved in the upstream folder.
% MODIFICATION HISTORY : 
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy, Seoul National University, 03/02/2020.
clear

[pth]=uigetfile_n_dir();  % Select the .tif file folder taken from the olympus Two photon microscope. 
ch=2;  % number of channel
ref=1; % Set the reference channel. For example, if you want to use the channel 1 as reference to calculate the motion artifact, 
       % NoRMCorr caculates the motion correction from ch1 and apply the shift to ch2.
step=30; % Number of frames in each stack. Usally resonant scanner scans in 30 fps. Therefore, 2 sec interval btw the stacks -> 60 frm
%%
for p=1:length(pth)
if ch~=1 && ref==2
    ap=1;
else
    ap=2;
end
clear order_img order
    img=imageDatastore(pth{1,p});
    sp=split(img.Files,'T');
    sp=split(sp(:,2),'.');
    for z=1:size(sp,1)
    order_img(z,1)=str2num(sp{z,1});
    end
    [as order]=sort(order_img,'ascend');
    %step=floor(size(im{1},3)/slice);
    slice=floor(size(order_img,1)/step/ch);
    for j=1:slice
    switch ch
        case 1
    im{1}=[]; im{2}=[];
    for z=1:step
        im{1}(:,:,z)= (char(img.Files(order((step*(j-1)+z),1),:)));
    end
    Y = single(im{1});                 % convert to single precision 
    T = size(Y,ndims(Y));
    Y = Y - min(Y(:));
    options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200);
    tic; [M1,shifts1,template1,options_rigid] = normcorre(Y,options_rigid); toc
    imwrite(uint16(mean(M1,3)),[pth{1,p}(1,1:end-7) '.tif'],'Writemode','Append')
        
        case 2
    im{1}=[]; im{2}=[];
    for z=1:step
        im{1}(:,:,z)=imread(char(img.Files(order(2*(step*(j-1)+z)-1,1),:)));
        im{2}(:,:,z)=imread(char(img.Files(order(2*(step*(j-1)+z),1),:)));
    end
    Y = single(im{ref});                 % convert to single precision 
    T = size(Y,ndims(Y));
    Y = Y - min(Y(:));
    options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',30,'us_fac',50,'init_batch',200);
    tic; [M1,shifts1,template1,options_rigid] = normcorre(Y,options_rigid); toc
    M2 = apply_shifts(im{ap},shifts1,options_rigid);
    I{1}=mean(M2,3); I{2}=mean(M1,3);
    s=subtraction_image_ftn_cell(I,0.7); %%%subtraction
    imwrite(uint16(mean(M2,3)),[pth{1,p}(1,1:end-7) '.tif'],'Writemode','Append')
    imwrite(uint16(mean(M1,3)),[pth{1,p}(1,1:end-7) '.tif'],'Writemode','Append')
    imwrite(uint16(s),[pth{1,p}(1,1:end-7) '_sub.tif'],'Writemode','Append')
    end
    end
end
%%

% if ischar(fnm)
%     fnm={fnm};
% end
% if ref==2
%     ap=1;
% else
%     ap=2;
% end
% for i=1:length(fnm)
%     nm=[fnm{1,i}(1,1:end-4) '_00.mat'];
%     if ~exist([pth nm]) % not exist
%         try
%     oir2stdData([pth fnm{1,i}],0,1);
%         catch
%     oir2stdData([pth fnm{1,i}],0,0);
%         end
%     end % exist
%         fnmss=dir(pth);
%         g=0;
%         for k=1:size(fnmss,1)
%             if contains(fnmss(k).name,fnm{1,i}(1,1:end-4)) && contains(fnmss(k).name,'.mat')
%             g=g+1;
%             arg(g,1)=k;
%             end
%         end
%     im{1}=[]; im{2}=[];
%     for k=1:size(arg,1)
%     load([pth fnmss(arg(k,1)).name])
%     im{1}(:,:,size(im{1},3)+1:size(im{1},3)+size(squeeze(stdData.Image{1}),3))=squeeze(stdData.Image{1});
%     im{2}(:,:,size(im{1},3)+1:size(im{1},3)+size(squeeze(stdData.Image{1}),3))=squeeze(stdData.Image{2});
%     end
%     %step=floor(size(im{1},3)/slice);
%     slice=floor(size(im{1},3)/step);
%     for j=1:slice
%     Y = single(im{ref}(:,:,step*(j-1)+1:step*j));                 % convert to single precision 
%     T = size(Y,ndims(Y));
%     Y = Y - min(Y(:));
%     options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200);
%     tic; [M1,shifts1,template1,options_rigid] = normcorre(Y,options_rigid); toc
%     M2 = apply_shifts(im{ap}(:,:,step*(j-1)+1:step*j),shifts1,options_rigid);
%     imwrite(uint16(mean(M2,3)),[pth fnm{1,i}(1,1:end-4) '.tif'],'Writemode','Append')
%     imwrite(uint16(mean(M1,3)),[pth fnm{1,i}(1,1:end-4) '.tif'],'Writemode','Append')
%     end
% end 