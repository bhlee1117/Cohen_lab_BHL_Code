function [f ax1]=show_Kymos_spikes(traces,otherT,cax)
noi=[1:size(traces,1)]; frmrate=1;
if nargin<4
    
    cax=[prctile(traces(:),1) prctile(traces(:),99)];
end
%figure;
%f=figure('units','normalized','outerposition',[0 0 1 1]);
%tiledlayout(11,1)
%ax1=nexttile([8 1]);
clf;
ax1=subplot(10,1,1:6);

t=[1:size(traces,2)]/frmrate; scale=max(prctile(traces,99,2))*3;
%tr=Result{i}.traces_res-median(Result{i}.traces_res,2); fprnt=Result{i}.c_ftprnt;
%tr=traces-movmedian(traces,150,2);
%tr=tr./get_threshold(tr,1);
%tr=(traces-median(traces,2,'omitnan'))./get_threshold(tr,1);
imagesc(traces,cax)
%S=ones(size(tr,1),size(tr,2)); S(~(spikes==1))=NaN;
% plot(t,S(noi,:)'+[1:size(noi,2)]*scale,'r.')
% set(gca,'ytick',[1:size(noi,2)]*scale,'yticklabel',noi)
%axis off
colormap(turbo)
ax2=subplot(10,1,7:10);

if size(otherT,1)>size(otherT,2)
    otherT=otherT';
end
for j=1:size(otherT,1)
plot(t,otherT(j,:))
hold all
end
%ylim([min(otherT(:)) max(otherT(:))])

linkaxes([ax1 ax2],'x')
end