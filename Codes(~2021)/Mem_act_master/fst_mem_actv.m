%%
clear

[ fnm pth]=uigetfile('*.tif','Select the TXN data','Multiselect','on');
[fnm_ca pth_ca ]=uigetfile('*.tif','Select the Calcium image');
fnm_ca={fnm_ca};
%%
    clear im imds subtraction imgs voldef_pl im_reg
imginf=imfinfo([pth_ca fnm_ca{1,1}]);
numstack=size(imginf,1);
for i=1:numstack/2
    im{1}(:,:,i)=imread([pth_ca fnm_ca{1,1}],2*i-1);
end

imgs = double(reshape(im{1}, [size(im{1}, 1), size(im{1}, 2), 1, 1, size(im{1},3)]))/4096;


opts = [];
opts.check_gradients = 155*0;
opts.mean_penalty = 1e-3;
opts.pix_resolution = [1,1];
opts.metric = 'nuclear';


opts.grid_spacing = [10, 10];
opts.isoTV = 5e-3;

opts.spline_order = 1;
opts.max_iters = 100;

opts.display = 'off';
opts.opt_method = 'lbfsg';

opts.border_mask = 6;
opts.k_down = 0.7;

tic
[voldef_pl, Tmin_pl,  Kmin_pl] = ptv_register(imgs, [], opts);
toc

im_reg=squeeze(voldef_pl);

savegif([char(pth_ca),char(fnm_ca),'.gif'], squeeze([imgs, voldef_pl]), 1/5);
for z=1:size(im_reg,3)
imwrite(uint16(im_reg(:,:,z)*4096), [pth_ca '\Reg_' char(fnm_ca) '.tif'],'tif','Writemode','Append','Compression', 'none')
end
%%
for p=1:size(fnm,2)
    clear im imds subtraction imgs voldef_pl im_reg
imginf=imfinfo([pth fnm{1,p}]);
numstack=size(imginf,1);
for i=1:numstack/2
    im{1}(:,:,i)=imread([pth fnm{1,p}],2*i-1);
    im{2}(:,:,i)=imread([pth fnm{1,p}],2*i);
end

subtraction=subtraction_image_ftn_cell(im);

imgs = double(reshape(subtraction, [size(subtraction, 1), size(subtraction, 2), 1, 1, size(subtraction,3)]))/4096;


opts = [];
opts.check_gradients = 155*0;
opts.mean_penalty = 1e-3;
opts.pix_resolution = [1,1];
opts.metric = 'nuclear';


opts.grid_spacing = [10, 10];
opts.isoTV = 5e-3;

opts.spline_order = 1;
opts.max_iters = 100;

opts.display = 'off';
opts.opt_method = 'lbfsg';

opts.border_mask = 6;
opts.k_down = 0.7;

tic
[voldef_pl, Tmin_pl,  Kmin_pl] = ptv_register(imgs, [], opts);
toc

im_reg=squeeze(voldef_pl);

%sp_pth=split(pth{1,p},'\');
savegif([pth,char(fnm{1,p}),'.gif'], squeeze([imgs, voldef_pl]), 1/5);
for z=1:size(im_reg,3)
imwrite(uint16(im_reg(:,:,z)*4096), [pth 'Sub_' char(fnm{1,p}) '.tif'],'tif','Writemode','Append','Compression', 'none')
end
end