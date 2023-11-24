function [f ax1]=show_traces_spikes(traces,spikes,otherT,frmrate)
if nargin<4
    frmrate=1;
end
%figure;
f=figure('units','normalized','outerposition',[0 0 1 1]);
%tiledlayout(11,1)
%ax1=nexttile([8 1]);
ax1=subplot(10,1,1:8);
noi=[1:size(traces,1)];
t=[1:size(traces,2)]/frmrate; scale=20;
%tr=Result{i}.traces_res-median(Result{i}.traces_res,2); fprnt=Result{i}.c_ftprnt;
tr=traces-movmedian(traces,150,2);
%tr=tr./get_threshold(tr,1);
tr=(traces-median(traces,2,'omitnan'))./get_threshold(tr,1);

lines=plot(t,tr(noi,:)'+[1:size(noi,2)]*scale); line_color=turbo(length(noi));
arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(line_color,2))
%arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(line_color,2))
hold all
S=tr; S(~(spikes==1))=NaN;
%S=ones(size(tr,1),size(tr,2)); S(~(spikes==1))=NaN;
plot(t,S(noi,:)'+[1:size(noi,2)]*scale,'r.')
set(gca,'ytick',[1:size(noi,2)]*scale,'yticklabel',noi)
%axis off

ax2=subplot(10,1,9:10);

if size(otherT,1)>size(otherT,2)
    otherT=otherT';
end
for j=1:size(otherT,1)
plot(t,otherT(j,:))
hold all
end
ylim([min(otherT(:)) max(otherT(:))])

linkaxes([ax1 ax2],'x')
end