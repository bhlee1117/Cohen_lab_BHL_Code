%%
clear
[fnm pth]=uigetfile('*.mat','Multiselect','On');

%%
ptls_total=[];
for i=1:length(fnm)
   
load([pth fnm{1,i}])
ptls=cell2mat(refinementData);
ptls_total=[ptls_total; ptls];
end

%%
pix=4;
storm_coor=ceil(ptls_total(:,1:2)*pix);
STORM_im=zeros(movieSize(1,1)*pix,movieSize(1,1)*pix);
for i=1:size(storm_coor,1)
STORM_im(storm_coor(i,1),storm_coor(i,2))=STORM_im(storm_coor(i,1),storm_coor(i,2))+1;
end

%%
imshow(STORM_im,[0 10])
colormap('gray')
axis equal tight off