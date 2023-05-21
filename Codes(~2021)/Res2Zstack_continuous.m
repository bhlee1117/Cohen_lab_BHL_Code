%%
%     Website : http://neurobiophysics.snu.ac.kr/
% This code makes averaged image from resonant scanning image. 
% Resonant scanned image -> Motion correction (NoRMCorr) -> average -> save in tiff 

% INPUTS 
% Olympus two photon images, oir format. 
% As oir2stdData library cannot import image more than 2 channel, import
% the image with 1 or 2 channel.
% 

% OUTPUTS
% one averaged image per one oir file
% MODIFICATION HISTORY : 
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy, Seoul National University, 03/02/2020.
clear

[fnm pth]=uigetfile('*.oir','Multiselect','On');  % Select the .oir file taken from the olympus Two photon microscope. 
ref=2; % Set the reference channel. For example, if you want to use the channel 1 as reference to calculate the motion artifact, 
       % NoRMCorr caculates the motion correction from ch1 and apply the shift to ch2.
slice=41;
gcp;
%%
if ref==2
    ap=1;
else
    ap=2;
end
for i=1:length(fnm)
    nm=[fnm{1,i}(1,1:end-4) '_00.mat'];
    if ~exist([pth nm]) % not exist
        try
    oir2stdData([pth fnm{1,i}],0,1);
        catch
    oir2stdData([pth fnm{1,i}],0,0);
        end
    else % exist
        fnmss=dir(pth);
        g=0;
        for k=1:size(fnmss,1)
            if contains(fnmss(k).name,fnm{1,i}(1,1:end-4)) && contains(fnmss(k).name,'.mat')
            g=g+1;
            arg(g,1)=k;
            end
        end
    end
    im{1}=[]; im{2}=[];
    for k=1:size(arg,1)
    load([pth fnmss(arg(k,1)).name])
    im{1}(:,:,size(im{1},3)+1:size(im{1},3)+size(squeeze(stdData.Image{1}),3))=squeeze(stdData.Image{1});
    im{2}(:,:,size(im{1},3)+1:size(im{1},3)+size(squeeze(stdData.Image{1}),3))=squeeze(stdData.Image{2});
    end
    
    for j=1:floor(size(im{1},3)/slice)
    Y = single(im{ref}(:,:,slice*(j-1)+1:slice*j));                 % convert to single precision 
    T = size(Y,ndims(Y));
    Y = Y - min(Y(:));
    options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200);
    tic; [M1,shifts1,template1,options_rigid] = normcorre(Y,options_rigid); toc
    M2 = apply_shifts(im{ap}(:,:,slice*(j-1)+1:slice*j),shifts1,options_rigid);
    imwrite(uint16(mean(M2,3)),[pth fnm{1,i}(1,1:end-4) '.tif'],'Writemode','Append')
    imwrite(uint16(mean(M1,3)),[pth fnm{1,i}(1,1:end-4) '.tif'],'Writemode','Append')
    end
end 