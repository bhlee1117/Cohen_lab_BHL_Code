clear;
[fnm pth]=uigetfile('*.tif','Select registerated Calcium image','Multiselect','on');
 %% load file
 if ischar(fnm)
     fnm={fnm};
 end
 for p=1:length(fnm)
im_inform=imfinfo([pth fnm{1,p}]);
for i=1:numel(im_inform)
    Y(:,:,i)=imread([pth fnm{1,p}],i);
end
Y = Y - min(Y(:)); 
if ~isa(Y,'single');    Y = single(Y);  end         % convert to single

[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;             
%total number of pixels
options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',30,'us_fac',50,'init_batch',200);
    tic; [Y,shifts1,template1,options_rigid] = normcorre(Y,options_rigid); toc
% Set parameters
pp=1;

[P,Y] = preprocess_data(Y,pp);
refined_data=manual_masking(Y);
 save([pth fnm{1,p}(1:end-4) '_astro_rst.mat'],'Y','refined_data','-v7.3');
 end