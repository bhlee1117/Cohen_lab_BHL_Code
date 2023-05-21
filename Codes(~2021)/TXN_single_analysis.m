clear
[fnm pth]=uigetfile('*.tif','Multiselect','on');
[fnm_dat pth_dat]=uigetfile('*.mat','Multiselect','on');
%%

for p=12%:length(fnm)
    clear im reslice TXN single
    th=[90 130;34 123];
    for j=1:2
inff=imfinfo([pth fnm{1,2*(p-1)+j}]);
load([pth_dat fnm_dat{1,2*(p-1)+j}])
for i=1:numel(inff)
    im(:,:,i)=imread([pth fnm{1,2*(p-1)+j}],i);
    reslice(:,i,:)=im(:,:,i);
end
threshold=th(j,:);
datum_filter=datum;
%datum_filter(find((datum(:,4)<100 & datum(:,9)>1000) | datum(:,3)<0),:)=NaN;
datum_filter(find((datum(:,4)<100 & datum(:,9)>1000) | datum(:,4)>1000 | datum(:,3)<0),:)=NaN;

TXN=find(datum_filter(:,4)>threshold(1,1) & datum_filter(:,7)>threshold(1,2));
single=find(datum_filter(:,4)<threshold(1,1) | datum_filter(:,7)<threshold(1,2));
CDSPP{p,j}=datum(TXN,[1 2 4 7]);

intensity=datum_filter(:,4).*datum_filter(:,7).*datum_filter(:,8).*datum_filter(:,9);
save(['H:\Image_data\KIST_RNAscope\40x\20190209_SNU_RNAscope\RNA_scope_confocal\Result\Coord_included_20200129\mRNA_' fnm{1,2*(p-1)+j} ,'.mat'],'threshold','single','TXN','intensity','datum')
RGB(:,:,j)=3*max(im,[],3);
    end
f = figure;
RGB(:,:,3)=0;
imagesc(RGB)
axis equal tight off
hold all
plot(CDSPP{p,1}(:,2)/160+1,CDSPP{p,1}(:,1)/160+1,'bo','markersize',15)
plot(CDSPP{p,2}(:,2)/160+1,CDSPP{p,2}(:,1)/160+1,'b.','markersize',15)
plot(datum(single,2)/160+1,datum(single,1)/160+1,'color',[0 1 1],'marker','.','markersize',10,'linestyle','none')

end