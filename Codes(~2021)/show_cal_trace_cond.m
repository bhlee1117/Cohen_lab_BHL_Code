%% This code is written to visualize calcium traces
% MODIFICATION HISTORY :
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy,
%           Seoul National University, 2020/08/14
%           Modified C_df, 2020/09/28
%           Add baseline correction from both starting part and end part,
%           2020/10/1
 function list=show_cal_trace_cond(dat,groups,m,post_ref,rep_n,cmap)

clear list
if post_ref
    Arc_class{m}=[dat{m}.Arc_post_ref(:,2:end)];
else
    Arc_class{m}=[dat{m}.Arc(:,2:end)];
end
for g=1:size(groups,2)
    for cd=1:size(groups{g},2)/2
        list{m,g}(:,cd)=Arc_class{m}(:,groups{g}(1,2*cd-1))==groups{g}(1,2*cd);
    end
    [rr cc]=cellfun(@size,dat{m}.Cal(:,1:3),'UniformOutput',false);
    list{m,g}=[list{m,g} sum(cell2mat(cc)>1,2)==3];
    list{m,g}=find(sum(list{m,g},2)==(size(groups{g},2)/2+1)); % cells that match the conditions
end
%%
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 1200 600]);
sc=1.5;
for d1=1:3
handles.axes1 = axes('Units','pixels','Position',[50+350*(d1-1) 50 300 500]);
gr=1*sc;
for g=1:size(groups,2)
     %cmap_blur=cmap(g,:)+0.2; cmap_blur(cmap_blur>1)=1;
    for n=rep_n{g}
        plot([1:size(dat{m}.run{d1},1)]/30,dat{m}.Cal{list{m,g}(n,1),d1}(1:size(dat{m}.run{d1},1))+gr,'color',cmap(g,:),'linewidth',1)
        hold all
%         for tr=1:size(dat{m}.Cal{list{m,g}(n,1),6+d1},1)
%         s=dat{m}.Cal{list{m,g}(n,1),6+d1}(tr,:);
%         if sum(s>size(dat{m}.run{d1},1))==0
%         plot([s(1,1):s(1,2)]/30,dat{m}.Cal{list{m,g}(n,1),d1}(1,s(1,1):s(1,2))+gr,'color',cmap(g,:),'linewidth',1)    
%         end
%         end
        gr=gr+1*sc;
    end
    %gr=gr+1*sc;
end
axis tight off
end
end