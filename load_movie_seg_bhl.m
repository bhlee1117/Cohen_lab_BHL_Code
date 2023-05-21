function m=load_movie_seg_bhl(mov_folder,dim,f_seg,target)

tic;
    for i=1:size(f_seg,2)-1
        t=ismember([f_seg(i):f_seg(i+1)-1],target);
        if sum(t)>0
          ind{i}=find(t);
        else
           ind{i}=[];
        end
    end

m=[];
for i=1:size(f_seg,2)-1
fname=['mc' num2str(i,'%02d') '.bin'];
if ~isempty(ind{i})
m_s=double(readBinMov_times(fullfile(mov_folder,fname),dim(1),dim(2),ind{i}));
m=cat(3,m,m_s);
end
end
toc;
end