%% Old population data: clusters and synchronization (Dragonfly rig, 2022.06 - 2023.02)
% Detailed population analysis on the results of Analysis_20261102_OldPopulationDataReanalyis.m
% (PopReanalysis_Result.mat next to each voltage movie: traces, spikes, complex spikes)
% 1. Load one recording and its analysis epoch (constant blue: blue on, onset transient left out)
% 2. Population clusters: hierarchical clustering of binned spike rates
% 3. Synchronization: co-activation events from spikes, CS onsets, subthreshold, trace
% Outputs (PopCluster_Result.mat, figures) are saved next to the voltage movie.
% Pilot: BHLpop120, constant blue, anesthetized, 142 cells (120 active), most high-SNR cells
% 2026.11.03, Byung Hun Lee

clear
clc;
cd '/Users/bhlee1117/Documents/GitHub/Cohen_lab_BHL_Code/Analysis_master_codes';
save_to='/Volumes/BHL18TB_D2/Arranged_Data/PopulationReanalysis';   % recList.mat of the re-analysis
recID='BHLpop120';     % pilot recording
figVis='on';           % 'off' to only save the figures

P=struct;
% analysis epoch
P.epoch='blue';        % 'blue': frames with blue on, 'all': whole recording
P.onsetMs=200;         % ms after each blue onset left out (onset transient)
P.minSpikes=5;         % cells with fewer spikes in the epoch are left out
% clustering
P.binMs=25; P.smoothBins=2;
P.maxClust=8;
P.minSil=0.2;          % mean silhouette below this -> one cluster (no structure)
% co-activation: synchronized population events, resolved by 4 signals
P.co.mods={'spike','CS','sub','trace'};
P.co.winMs=[10 30 30 10];      % coincidence window per signal
P.co.minSepMs=[30 60 60 30];   % min. interval between events
P.co.subThr=2.5;               % z, subthreshold depolarization counted as active
P.co.trThr=3.5;                % z, trace counted as active
P.co.minCells=3;               % min. co-active cells for an event
P.co.nShuf=20; P.co.shufPrc=99.5;   % event threshold: percentile of circularly shifted surrogates
P.co.nNull=200;                % null for population consistency (same event sizes, cell rates)
P.co.motifCut=0.6;             % 1-cosine cutoff grouping events into motifs
P.co.minMotif=3;               % events for a recurring motif
P.co.subLpMs=10; P.co.subBaseMs=300;   % subthreshold

%% 1. Load the recording and its analysis epoch
load(fullfile(save_to,'recList.mat'))
f=find(strcmp(recList.RecID,recID));
load(fullfile(recList.Path{f},'PopReanalysis_Result.mat'))
ep=analysis_epoch(Result,P);
act=find(sum(Result.spike(:,ep),2)>=P.minSpikes);
fprintf('%s %s | %s: %s | %s | epoch %.1f of %.1f s | %d/%d cells with >= %d spikes\n', ...
    Result.RecID,Result.recName,Result.protocol,Result.protocolDetail,Result.state, ...
    sum(ep)/Result.frm_rate,Result.nFrames/Result.frm_rate,numel(act),size(Result.spike,1),P.minSpikes);

%% 2. Population clusters
Clu=population_cluster(Result,act,ep,P);
plot_clusters(Result,Clu,ep,P,figVis);
fprintf('%d clusters (best mean silhouette %.2f)\n',Clu.k,max([Clu.sil -1]));

%% 3. Synchronization: co-activation events
rng(0,'twister');   % surrogates and null samples: the same result on every run
Co=coactivation(Result,act,ep,P);
plot_coactivation(Result,Co,Clu,P,figVis);
for m=1:numel(P.co.mods)
    E=Co.(P.co.mods{m});
    fprintf('%-6s %4d events (%.1f /s), consistency z %.1f, %d recurring motifs\n',P.co.mods{m}, ...
        E.nEvents,E.nEvents/(sum(ep)/Result.frm_rate),E.zCons,E.nMotif);
end
save(out_file(Result.fpath,'Result.mat'),'Clu','Co','ep','act','P','recID');

%% Functions
function fn=out_file(fpath,suffix)
% per-recording outputs are saved next to the voltage movie
fn=fullfile(fpath,['PopCluster_' suffix]);
end

function ep=analysis_epoch(Result,P)
% frames analyzed: not excluded (laser off, motion transients, global artifacts); with
% P.epoch='blue' only while blue is on, without the first P.onsetMs after each blue onset
ep=~Result.excluded(:)';
if strcmp(P.epoch,'blue') && any(Result.Blue>0)
    on=Result.Blue(:)'>0;
    nOn=round(P.onsetMs*1e-3*Result.frm_rate);
    for k=find(diff([false on])==1), on(k:min(end,k+nOn-1))=false; end
    ep=ep & on;
end
end

function C=population_cluster(Result,act,ep,P)
% hierarchical clustering (average linkage, 1-correlation) of binned spike rates in the
% epoch, number of clusters: smallest k within 90% of the best mean silhouette
nCell=size(Result.spike,1);
C=struct('cells',act,'label',zeros(nCell,1),'order',act,'R',[],'k',0,'sil',[], ...
    'rate',[],'tBin',[],'okBin',[],'Z',[]);
if numel(act)<4
    C.label(act)=1; C.k=double(~isempty(act)); return
end
binF=max(1,round(P.binMs*1e-3*Result.frm_rate));
nB=floor(Result.nFrames/binF);
rate=squeeze(sum(reshape(double(Result.spike(act,1:nB*binF)),numel(act),binF,nB),2));
okB=all(reshape(ep(1:nB*binF),binF,nB),1);
rate=smoothdata(rate,2,'gaussian',P.smoothBins*5)*Result.frm_rate/binF;
R=corr(rate(:,okB)'); R(~isfinite(R))=0;
D=1-R; D=(D+D')/2; D(1:numel(act)+1:end)=0;
Dv=squareform(D,'tovector');
Z=linkage(Dv,'average');
kMax=min(P.maxClust,numel(act)-1);
sil=nan(1,kMax);
for k=2:kMax
    lab=cluster(Z,'maxclust',k);
    s=silhouette(rate(:,okB),lab,'correlation');
    cnt=accumarray(lab,1);
    s(cnt(lab)==1)=0;        % singletons count as 0 (MATLAB returns 1)
    sil(k)=mean(s,'omitnan');
end
if max(sil)<P.minSil
    k=1; lab=ones(numel(act),1);
else
    k=find(sil>=0.9*max(sil),1);
    lab=cluster(Z,'maxclust',k);
end
leaf=optimalleaforder(Z,Dv);
[~,bySize]=sort(accumarray(lab,1),'descend');   % cluster 1 = largest
newLab=zeros(size(lab));
for q=1:k, newLab(lab==bySize(q))=q; end
pos=zeros(numel(act),1); pos(leaf)=1:numel(act);
[~,o]=sortrows([newLab pos]);
C.label(act)=newLab; C.order=act(o); C.R=R(o,o); C.k=k; C.sil=sil; C.Z=Z;
C.rate=rate(o,:); C.tBin=((1:nB)-0.5)*binF/Result.frm_rate; C.okBin=okB;
end

function co=coactivation(Result,act,ep,P)
% synchronized population events among the active cells in the epoch, from 4 signals:
%   spike, CS (complex-spike onsets), sub (subthreshold, spikes removed), trace (z-trace)
% per event: which cells take part; across events: does the same population recur
f=@(ms) max(1,round(ms*1e-3*Result.frm_rate));
valid=ep & ~imdilate([false diff(Result.Blue(:)'>0)~=0],ones(1,2*f(4)+1));  % blue on/off edges
A=struct;
A.spike=Result.spike(act,:);
A.CS=Result.CSonset(act,:);
tr=double(Result.dFF(act,:)); tr(~isfinite(tr))=0;
tr(imdilate(Result.spike(act,:),ones(1,2*f(3)+1)))=NaN;   % spikes +-3 ms removed
tr=fillmissing(tr,'linear',2,'EndValues','nearest');
sub=movmean(tr,f(P.co.subLpMs),2); sub=sub-movmedian(sub,f(P.co.subBaseMs),2);
sub=sub./(1.4826*mad(sub(:,valid),1,2)+eps);
A.sub=sub>P.co.subThr;
A.trace=movmean(double(Result.zTr(act,:)),f(2),2)>P.co.trThr;
for m=1:numel(P.co.mods)
    md=P.co.mods{m};
    E=coact_events(A.(md),valid,f(P.co.winMs(m)/2),f(P.co.minSepMs(m)),P);
    E.cells=act;
    if E.nEvents>0, E=coact_consistency(E,P); end
    co.(md)=E;
end
end

function E=coact_events(A,valid,w,sep,P)
% event = peak of the number of cells active within +-w frames, above the
% shufPrc percentile of circularly shifted surrogates (and >= minCells)
E=struct('nEvents',0,'frame',[],'nCells',[],'thr',NaN,'part',false(size(A,1),0), ...
    'motif',[],'nMotif',0,'motifFrac',NaN,'motifPart',[], ...
    'consistency',NaN,'nullMean',NaN,'zCons',NaN,'pCons',NaN);
[nC,nT]=size(A);
if nC<P.co.minCells || ~any(A(:)), return; end
A(:,~valid)=false;
D=movmax(double(A),2*w+1,2)>0;
c=sum(D,1);
cs=zeros(P.co.nShuf,nT,'single');
for s=1:P.co.nShuf
    sh=randi(nT,nC,1); Ds=false(nC,nT);
    for n=1:nC, Ds(n,:)=circshift(D(n,:),sh(n)); end
    cs(s,:)=sum(Ds,1);
end
E.thr=max(P.co.minCells,floor(prctile(cs(:),P.co.shufPrc))+1);
if max(c)<E.thr, return; end
[~,loc]=findpeaks(double(c),'MinPeakHeight',E.thr-0.5,'MinPeakDistance',sep);
E.nEvents=numel(loc); E.frame=loc;
E.part=D(:,loc); E.nCells=sum(E.part,1);
end

function E=coact_consistency(E,P)
% mean pairwise cosine similarity of participation vectors vs. a null with the
% same event sizes and per-cell participation rates; events grouped into motifs
X=double(E.part); [nC,nE]=size(X);
E.motif=ones(1,nE);
if nE<2, return; end
sz=sum(X,1);
mcos=@(Y) (sum(sum(Y./sqrt(sum(Y,1)),2).^2)-nE)/(nE*(nE-1));   % binary columns
E.consistency=mcos(X);
rate=max(mean(X,2),1e-6);
nul=zeros(P.co.nNull,1);
for s=1:P.co.nNull   % weighted sampling without replacement (Efraimidis-Spirakis keys)
    [~,ord]=sort(rand(nC,nE).^(1./rate),1,'descend');
    Rk=zeros(nC,nE); Rk(sub2ind([nC nE],ord,repmat(1:nE,nC,1)))=repmat((1:nC)',1,nE);
    nul(s)=mcos(double(Rk<=sz));
end
E.nullMean=mean(nul); E.zCons=(E.consistency-mean(nul))/(std(nul)+eps);
E.pCons=(sum(nul>=E.consistency)+1)/(P.co.nNull+1);

idx=1:nE;
if nE>5000, idx=sort(randperm(nE,5000)); end
lab=cluster(linkage(X(:,idx)','average','cosine'),'cutoff',P.co.motifCut,'criterion','distance');
if nE>5000   % remaining events -> nearest motif centroid
    Cm=zeros(nC,max(lab)); for q=1:max(lab), Cm(:,q)=mean(X(:,idx(lab==q)),2); end
    [~,labAll]=max((X./sqrt(sum(X,1)))'*(Cm./sqrt(sum(Cm.^2,1))),[],2);
    labAll(idx)=lab; lab=labAll;
end
cnt=accumarray(lab(:),1);
[~,bySize]=sort(cnt,'descend');
newLab=zeros(1,nE);
for q=1:numel(bySize), newLab(lab==bySize(q))=q; end
E.motif=newLab;
E.nMotif=sum(cnt>=P.co.minMotif);
E.motifFrac=mean(ismember(newLab,1:E.nMotif));
E.motifPart=zeros(nC,min(E.nMotif,10));
for q=1:size(E.motifPart,2), E.motifPart(:,q)=mean(X(:,newLab==q),2); end
end

function plot_clusters(Result,C,ep,P,figVis)
% footprints and spike raster colored by cluster, rate correlation, silhouette per k
nCell=size(Result.ftprnt,3);
cmapC=lines(max(C.k,1));
col=repmat([0.75 0.75 0.75],nCell,1);
col(C.label>0,:)=cmapC(C.label(C.label>0),:);
t=(1:Result.nFrames)/Result.frm_rate;

fig=figure(310); clf; theme(fig,'light');
set(fig,'Color','w','Position',[50 50 1600 1100],'Visible',figVis);
tiledlayout(10,6,'TileSpacing','compact','Padding','compact');
nexttile(1,[4 2]);
fp=Result.ftprnt>0.3*max(Result.ftprnt,[],[1 2]);
show_footprnt_contour(fp,Result.ref_im,col);
title(sprintf('%d footprints, %d in clustering, %d clusters',nCell,numel(C.cells),C.k));

nexttile(3,[4 2]);
if ~isempty(C.R)
    imagesc(C.R,[-0.2 0.8]); axis image; colormap(gca,turbo); colorbar;
    hold on
    b=cumsum(accumarray(C.label(C.order),1));
    for q=1:numel(b)-1, xline(b(q)+0.5,'w'); yline(b(q)+0.5,'w'); end
    set(gca,'xtick',[],'ytick',[]);
    title(sprintf('Spike-rate correlation (%d ms bins, sorted by cluster)',P.binMs));
else
    text(0.1,0.5,'< 4 cells, no clustering'); axis off
end

nexttile(5,[2 2]);
if ~isempty(C.sil)
    plot(2:numel(C.sil),C.sil(2:end),'ko-','MarkerFaceColor','k'); hold on
    yline(P.minSil,'r--'); xline(C.k,'b:');
    xlabel('# clusters'); ylabel('mean silhouette');
    title(sprintf('k = %d (below %.1f: one cluster)',C.k,P.minSil),'FontSize',9);
else
    axis off
end
nexttile(17,[2 2]); hold on
if C.k>0 && ~isempty(C.rate)
    r=C.rate; r(:,~C.okBin)=NaN;
    for q=1:C.k
        m=mean(r(C.label(C.order)==q,:),1);
        plot(C.tBin,m/max(m+eps)+q-1,'Color',cmapC(q,:));
    end
    axis tight; xlabel('time (s)'); ylabel('cluster'); title('Cluster mean rate (norm.)','FontSize',9);
else
    axis off
end

ax1=nexttile(25,[1 6]);
binF=max(1,round(10e-3*Result.frm_rate));
nB=floor(Result.nFrames/binF);
pr=sum(reshape(double(Result.spike(C.cells,1:nB*binF)),numel(C.cells),binF,nB),[1 2]);
bar((1:nB)*binF/Result.frm_rate,squeeze(pr),1,'k','EdgeColor','none'); hold on
yl=ylim; area(t,~ep*yl(2),'FaceColor',[0.85 0.85 0.85],'EdgeColor','none','FaceAlpha',0.5);
ylabel('spikes / 10 ms'); axis tight; set(gca,'xticklabel',[]);
title(sprintf('%s  %s  |  %s: %s  |  %s  |  gray: outside the analysis epoch', ...
    Result.RecID,Result.recName,Result.protocol,Result.protocolDetail,Result.state), ...
    'Interpreter','none','FontSize',9);

ax2=nexttile(31,[5 6]); hold on
for k=1:numel(C.order)
    n=C.order(k);
    s=find(Result.spike(n,:) & ~Result.CStrace(n,:));
    plot(t(s),k*ones(size(s)),'.','Color',col(n,:),'MarkerSize',5);
    s=find(Result.spike(n,:) & Result.CStrace(n,:));
    plot(t(s),k*ones(size(s)),'r.','MarkerSize',5);
end
set(gca,'YDir','reverse'); axis tight; xlabel('time (s)'); ylabel('cells (cluster order)');
title('Spikes colored by cluster, red: spikes in complex spikes','FontSize',9);
linkaxes([ax1 ax2],'x');
exportgraphics(fig,out_file(Result.fpath,'clusters.png'),'Resolution',150);
end

function plot_coactivation(Result,co,C,P,figVis)
mods=P.co.mods;
names={'Spike','Complex spike','Subthreshold','Trace'};
t=(1:Result.nFrames)/Result.frm_rate;
fig=figure(311); clf; theme(fig,'light');
set(fig,'Color','w','Position',[50 50 1800 1250],'Visible',figVis);
tiledlayout(4,4,'TileSpacing','compact','Padding','compact');
for m=1:4
    E=co.(mods{m});
    if E.nEvents==0
        nexttile(m); text(0.05,0.5,sprintf('%s: no synchronized events',names{m})); axis off
        for r=2:4, nexttile((r-1)*4+m); axis off; end
        continue
    end
    nM=max(E.nMotif,1);
    cm=[lines(nM); 0.6 0.6 0.6];
    ec=min(E.motif,nM+1); ec(E.motif>E.nMotif)=nM+1;

    nexttile(m); hold on
    scatter(t(E.frame),E.nCells,18,cm(ec,:),'filled');
    yline(E.thr,'k--');
    xlim(t([1 end])); xlabel('time (s)'); ylabel('# co-active cells');
    title(sprintf('%s: %d events',names{m},E.nEvents));

    show=1:E.nEvents;
    if E.nEvents>400, show=sort(randperm(E.nEvents,400)); end
    [~,o]=sortrows([E.motif(show)' E.frame(show)']); show=show(o);
    X=double(E.part(:,show)); Xn=X./sqrt(sum(X,1));
    nexttile(4+m);
    imagesc(Xn'*Xn,[0 1]); axis image; colormap(gca,turbo); colorbar;
    title(sprintf('event similarity (cosine)  mean %.2f vs null %.2f, z=%.1f',E.consistency, ...
        E.nullMean,E.zCons),'FontSize',8);

    nexttile(8+m);
    [~,cellOrd]=ismember(C.order(:),E.cells); cellOrd=cellOrd(cellOrd>0);
    imagesc(X(cellOrd,:)); colormap(gca,flipud(gray));
    hold on; b=find(diff(E.motif(show))~=0);
    for q=b, xline(q+0.5,'r'); end
    xlabel('events (sorted by motif)'); ylabel('cells (cluster order)');
    title(sprintf('%d recurring motifs (>=%d events) hold %d%% of events',E.nMotif, ...
        P.co.minMotif,round(100*E.motifFrac)),'FontSize',8);

    nexttile(12+m);
    imshow2(Result.ref_im,[]); hold on
    xy=Result.centroid(E.cells,:);
    scatter(xy(:,1),xy(:,2),10,[1 1 1]*0.8);
    for q=1:min(3,size(E.motifPart,2))
        pq=E.motifPart(:,q); s=pq>0.3;
        scatter(xy(s,1),xy(s,2),20+80*pq(s),cm(q,:),'filled','MarkerFaceAlpha',0.8);
    end
    title('cells in motif 1-3 (size = participation)','FontSize',8);
end
sgtitle(sprintf('%s  %s  |  %s: %s  |  co-activation of %d cells in the analysis epoch',Result.RecID, ...
    Result.recName,Result.protocol,Result.protocolDetail,numel(C.cells)), ...
    'Interpreter','none','FontSize',10);
exportgraphics(fig,out_file(Result.fpath,'coactivation.png'),'Resolution',120);
end
