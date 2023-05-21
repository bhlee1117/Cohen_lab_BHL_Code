%%
clear

pth=uigetfile_n_dir();
%%
for p=1:size(pth,2)
    clear im imds subtraction imgs voldef_pl im_reg
imginf=imfinfo(pth{1,p});
numstack=size(imginf,1);
for i=1:numstack/2
    im{1}(:,:,i)=imread(pth{1,p},2*i-1);
    im{2}(:,:,i)=imread(pth{1,p},2*i);
end

subtraction=subtraction_image_ftn_cell_dend(im,1.1);

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

sp_pth=split(pth{1,p},'\');
savegif(['H:\Image_data\OlympusTPM\20191129_OE51_CtxAA_dend\Ptv_reg\',char(sp_pth(end,1)),'.gif'], squeeze([imgs, voldef_pl]), 1/5);
for z=1:size(im_reg,3)
imwrite(uint16(im_reg(:,:,z)*4096), ['H:\Image_data\OlympusTPM\20191205_jRCaMP1a\jRCAmp' '\Sub_' char(sp_pth(end,1)) '.tif'],'tif','Writemode','Append','Compression', 'none')
end
end