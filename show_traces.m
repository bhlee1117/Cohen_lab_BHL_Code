function show_traces(traces_hi,scale,spike,Blue)


traces_hi=traces_hi-median(traces_hi,2);
n=size(traces_hi,1);
if nargin<3
    spike=[];
end

if nargin<4
    Blue=[];
    figure;
tiledlayout(7,4)
ax2 = nexttile([7 4]);
lines=plot(traces_hi'+[1:size(traces_hi,1)]*scale);
arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(jet(n),2))
axis tight
set(gca,'ytick',[1:size(traces_hi,1)]*scale,'yticklabel',[1:size(traces_hi,1)])
S=traces_hi; S(~spike)=NaN;
plot(S'+[1:size(traces_hi,1)]*scale,'r.')

else

figure;
tiledlayout(10,4)
% ax1 = nexttile([2 2]);
% colr = max(colormap(jet(size(c_ftprnt,3))),0);
% imshow2(squeeze(sum(c_ftprnt.*reshape(colr,1,1,[],3),3)),[]); hold all;
% text(coord(:,1)',coord(:,2)',num2str([1:size(c_ftprnt,3)]'),'color','w')
% 
% ax4 = nexttile([2 2]);
% imshow2(mean(mov_mc,3),[]); colormap('gray')

%linkaxes([ax1 ax4],'xy')

ax2 = nexttile([7 4]);
lines=plot(traces_hi'+[1:size(traces_hi,1)]*scale);
arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(jet(n),2))
axis tight
set(gca,'ytick',[1:size(traces_hi,1)]*scale,'yticklabel',[1:size(traces_hi,1)])
S=traces_hi; S(~spike)=NaN;
plot(S'+[1:size(traces_hi,1)]*scale,'r.')

ax3 = nexttile([1 4]);
plot(Blue)
end