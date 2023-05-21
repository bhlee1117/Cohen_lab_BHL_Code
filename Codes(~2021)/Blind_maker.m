clear
[fnm pth]= uigetfile('.tif','Select the images to be blindly transformed','multiselect','on');

%%
ss=randn(1,size(fnm,2));
[a answer]=sort(ss,'descend');
mkdir([pth 'Blind'])
for i=1:size(fnm,2)
  sp=split(fnm{1,i},'_');
  copyfile([pth fnm{1,i}],[pth 'Blind' '\' sp{1,1} '_' num2str(answer(1,i)) '.tif'])
  Blind_answer{1,i}=sp{2,1};
Blind_answer{2,i}=answer(1,i);  
end
%%
save([pth '\Blind_answer.mat'],'Blind_answer')
