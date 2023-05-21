
m=2; day=1;
[r c]=cellfun(@size,dat{m}.Cal(:,1));
rng=find(dat{m}.Arc_post_ref(:,2)==2 & dat{m}.Arc_post_ref(:,3)==2 & c>1);
rng([6 14])=[];
regist=Full_result{1, m}.regist;
options=Full_result{m}.options;
cmap=[1 0.2 0];
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 500 500]);
for t=1:size(Y,3)
clf('reset')    
Y_crop=imresize(imcrop(imrotate(imresize(Y(:,:,t),2),regist.angle(day,1),'crop'),regist.roi(day,:)*2),0.5);
imagesc(Y_crop)
axis equal tight off
hold all
colormap('gray')
for i=rng'
A_temp = full(reshape(Full_result{m}.SP_identified{i,day},options.d1,options.d2));
A_temp = medfilt2(A_temp,[3,3]);
A_temp = A_temp(:);
 [temp,ind] = sort(A_temp(:).^2,'ascend'); 
 temp =  cumsum(temp);
 ff = find(temp > (1-options.nrgthr)*temp(end),1,'first');
 p=full(reshape(Full_result{m}.SP_identified{i,day},options.d1,options.d2));
 p(p > (1-options.nrgthr)*temp(end))=0;
CC{i} = contour(reshape(A_temp,options.d1,options.d2),[0,0]+A_temp(ind(ff)),'LineColor',cmap(1,:), 'linewidth', 1);
end
F=getframe;
imwrite(uint8(F.cdata),['D:\Byunghun_Lee\VR\Images\Movie\Mov' num2str(t) '.jpeg']);
end