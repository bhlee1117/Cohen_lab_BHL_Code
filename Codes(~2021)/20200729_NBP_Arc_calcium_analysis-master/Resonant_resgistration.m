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
%           Update by Byung Hun Lee, Deptartment of Physics and Astronomy, Seoul National University, 08/12/2020.
clear

[pth]=uigetfile_n_dir();  % Select the .tif file folder taken from the olympus Two photon microscope. 
ch=3;  % number of channel
ref=2; % Set the reference channel. For example, if you want to use the channel 1 as reference to calculate the motion artifact, 
       % NoRMCorr caculates the motion correction from ch1 and apply the shift to ch2.
step=120; % Number of frames in each stack. Usally resonant scanner scans in 30 fps. Therefore, 2 sec interval btw the stacks -> 60 frm
zstack=true ;
%%
    f = waitbar(0,['Registration on progress']);
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
    if isempty(step)
    step=size(sp,1)/ch;
    end
    slice=floor(size(order_img,1)/step/ch);

    for j=1:slice
        waitbar((slice*(p-1)+j)/slice/length(pth),f,[num2str(p) '/' num2str(length(pth)) 'th image.' 'Registration on progress,' num2str(j) '/' num2str(slice) 'th slice.']);
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
    tic; [Y,shifts1,template1,options_rigid] = normcorre(Y,options_rigid); toc
    M1(:,:,j)=mean(Y,3);
    M2(:,:,j) =mean(apply_shifts(im{ap},shifts1,options_rigid),3);
    I{1}=M2(:,:,j); I{2}=M1(:,:,j);
    s(:,:,j)=subtraction_image_ftn_cell(I,0.9); %%%subtraction
%     imwrite(uint16(mean(M2,3)),[pth{1,p}(1,1:end-7) '.tif'],'Writemode','Append')
%     imwrite(uint16(mean(M1,3)),[pth{1,p}(1,1:end-7) '.tif'],'Writemode','Append')
%     imwrite(uint16(s),[pth{1,p}(1,1:end-7) '_sub.tif'],'Writemode','Append')
    end
    end
    
    if zstack
    else
    tic; [M1,shifts2,template1,options_rigid] = normcorre(M1,options_rigid); toc
    M2 = apply_shifts(M2,shifts2,options_rigid);
    s=apply_shifts(s,shifts2,options_rigid);
    end
    for j=1:slice
    imwrite(uint16(M1(:,:,j)),[pth{1,p}(1,1:end-7) '.tif'],'Writemode','Append')
    imwrite(uint16(M2(:,:,j)),[pth{1,p}(1,1:end-7) '.tif'],'Writemode','Append')
    imwrite(uint16(s(:,:,j)),[pth{1,p}(1,1:end-7) '_sub.tif'],'Writemode','Append')
    end
end
close(f)