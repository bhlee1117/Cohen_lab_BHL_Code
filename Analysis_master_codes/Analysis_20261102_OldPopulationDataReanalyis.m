%% Old population data re-analysis (Dragonfly rig, 2022.06 - 2023.02)
% 1. Expand the parent folders listed in PopulationData_Arrangement.xlsx (Sheet1)
%    into one row per recording, classify the protocol, flag recordings without
%    a motion-corrected movie (mc.bin) -> sheet 'PopulationRecordings'
% 2. Re-extract traces from the artifact-corrected movie
%    (mc.bin -> bleaching correction -> motion regression) with extractVoltronST
% 3. Spike and complex-spike detection, figures per recording
% 4. Summary across experiments, written back to the sheet
% Per-recording outputs (PopReanalysis_Result.mat, figures, movie) are saved next to the
% voltage movie; save_to keeps only recList.mat and the cross-recording summary.
% Clusters and synchronization: Analysis_20261103_OldPopulationCluster.m
% 2026.11.02, Byung Hun Lee

clear
clc;
cd '/Users/bhlee1117/Documents/GitHub/Cohen_lab_BHL_Code/Analysis_master_codes';
xlsFile=['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'PopulationData_Arrangement.xlsx'];
recSheet='PopulationRecordings';
save_to='/Volumes/BHL18TB_D2/Arranged_Data/PopulationReanalysis';   % recList + summary
if ~exist(save_to,'dir'), mkdir(save_to); end

interactiveUI=false;   % true: tune extractVoltronST with its sliders on the first detection window
csUI=false;            % true: tune complex-spike settings in interactive_cs_finder (most active cell)
figVis='on';           % 'off' to only save the figures
recToRun=[];           % rows of recList to process, [] = all usable recordings
overwrite=false;       % false: skip recordings that already have a result file

P=struct;
% artifact correction / extraction
P.blockLen=5000;       % frames per motion-regression / projection block
P.detWinLen=3000;      % frames per detection window
P.maxDetWin=3;         % detection windows, short recordings
P.maxDetWinLong=6;     % detection windows, recordings > 20000 frames
P.mergeRadius=4;       % px, footprints from different windows closer than this are merged
P.edgeMask=12;         % px, border band zeroed before detection (motion-correction edges)
P.detSmooth=3;         % frames, movmean of the detection movie (spikes span ~2-3 frames)
P.vopts=struct('hpWin',30,'smoothSigma',1,'zSeedMin',3.5,'zThrFoot',2.5, ...
    'minPeakDist_fr',10,'minPeakProm',0.75,'minArea',12,'seedStride',2, ...
    'r_max',6,'bound',12,'minEventsSyn',3,'footSigma',1,'usePeakWeight',true);
P.maxArea=400;         % px, larger footprints are dropped (vessels, neuropil)
P.artZ=6;              % robust z of frame-to-frame RMS change -> fast-motion artifact frame
P.artDilate=10;        % frames excluded around each artifact frame
P.negZ=4;              % motion-like frame: many cells deflect below -negZ together
P.globalFrac=0.1;      % ... fraction of cells above max(this, median+6*MAD)
P.globalDilate=15;     % frames excluded around each such frame
% spikes
P.hiWin=100;           % frames, movmedian high-pass before spike detection
P.spThr=4.5; P.spProm=2;
P.minSpikes=5;         % cells with fewer spikes are not active
P.placeBin=50;
% complex spikes: interactive_cs_finder settings (envelope on/off threshold, # spikes,
% min. AUC in z*frames, # spikes for the mean ISI; ISI < 20 frames is fixed in the function)
P.cs=struct('thres',[4 1],'nSpike',3,'maxAUC',50,'nSpikeISI',3);
% spikes on a plateau: the 100-frame high-pass of the regular detection removes the plateau and
% stunts the spikes on it, so inside putative CS (envelope transients), after the first spike
% (the CS onset), peaks of z reaching P.spThr above the baseline in front of the plateau
% (median over baseFr before the first spike) with a 7-pt amplitude (the prominence of
% find_spike_bh) >= prom are added. Plateaus >= noisyRatio x noisier than the cell's quiet
% baseline (frame-to-frame differences over noiseWinFr, regular spikes left out): nothing is
% added, and regular spikes are kept only if their 7-pt amplitude is >= P.spProm x that ratio
P.cs.baseFr=15; P.cs.minIsiFr=3; P.cs.prom=3; P.cs.noiseWinFr=31; P.cs.noisyRatio=3;

%% 1. Recording list: one row per recording
raw=readcell(xlsFile,'Sheet','Sheet1');
isParent=cellfun(@(x) isa(x,'missing'),raw(:,1)) & ...
    cellfun(@(x) ischar(x) && startsWith(x,'/Volumes'),raw(:,3));
parentDirs=raw(isParent,3);

recList=[];
for p=1:numel(parentDirs)
    recDirs=find_recordings(parentDirs{p});
    for r=1:numel(recDirs)
        recList=[recList; inspect_recording(recDirs{r},parentDirs{p})];
    end
    fprintf('%s: %d recordings\n',parentDirs{p},numel(recDirs));
end
recList=struct2table(recList);
[~,o]=sortrows([recList.Date str2double(extractBefore(recList.Recording,7))]);
recList=recList(o,:);
recList.RecID=compose('BHLpop%03d',(1:height(recList))');
listFile=fullfile(save_to,'recList.mat');
if exist(listFile,'file')   % keep extraction times and summary columns of the previous list
    old=load(listFile); old=old.recList;
    [tf,io]=ismember(recList.Path,old.Path);
    for v=setdiff(old.Properties.VariableNames,recList.Properties.VariableNames,'stable')
        x=old.(v{1});
        if iscell(x), y=repmat({''},height(recList),1); else, y=nan(height(recList),1); end
        y(tf)=x(io(tf)); recList.(v{1})=y;
    end
end

backupFile=strrep(xlsFile,'.xlsx','_backup_20261102.xlsx');
if ~exist(backupFile,'file'), copyfile(xlsFile,backupFile); end
write_rec_sheet(recList,xlsFile,recSheet);
save(listFile,'recList');

disp(groupsummary(recList,{'Protocol','Use'}));
fprintf('%d recordings, %d usable, %d without mc.bin\n',height(recList), ...
    sum(recList.Use),sum(strcmp(recList.MCstatus,'NO MC')));

%% 2. Re-extract traces from the artifact-corrected movie
load(fullfile(save_to,'recList.mat'))
if isempty(recToRun), recToRun=find(recList.Use==1)'; end

for f=recToRun
    resFile=out_file(recList.Path{f},'Result.mat');
    if exist(resFile,'file') && ~overwrite, continue; end
    fprintf('Processing %s %s (%d/%d)...\n',recList.RecID{f},recList.Path{f}, ...
        find(recToRun==f),numel(recToRun));
    tic;
    try
        Result=reextract_recording(recList(f,:),P,interactiveUI);
        safe_save(resFile,Result);
        fprintf('  %d footprints, %.1f min\n',size(Result.ftprnt,3),toc/60);
    catch ME
        fprintf(2,'  FAILED %s: %s (%s line %d)\n',recList.RecID{f},ME.message, ...
            ME.stack(1).name,ME.stack(1).line);
    end
end

%% 3. Spikes, complex spikes, figures
for f=recToRun
    resFile=out_file(recList.Path{f},'Result.mat');
    if ~exist(resFile,'file'), continue; end
    try
        load(resFile)
    catch ME
        fprintf(2,'  UNREADABLE %s: %s (re-run section 2 for it)\n',recList.RecID{f},ME.message);
        continue
    end
    if isempty(Result.ftprnt), continue; end
    if ~isfield(Result,'extractedOn') && ismember('Extracted',recList.Properties.VariableNames)
        Result.extractedOn=recList.Extracted{f};   % results extracted before the time stamp existed
    end
    try
        Result=refresh_stim(Result);
        Result=analyze_recording(Result,P,csUI);
        safe_save(resFile,Result);
        plot_recording(Result,P,figVis);
        plot_cs_check(Result,P,figVis);
        fprintf('%s: %d/%d active cells, median %.1f Hz, %d CS (re-detected on plateaus: %d added, %d removed)\n', ...
            Result.RecID,sum(Result.isActive),numel(Result.isActive), ...
            median(Result.spRate(Result.isActive)),sum(Result.nCS),nnz(Result.spikeAdded),nnz(Result.spikeRemoved));
    catch ME
        fprintf(2,'  FAILED %s: %s (%s line %d)\n',Result.RecID,ME.message, ...
            ME.stack(1).name,ME.stack(1).line);
    end
end

%% 4. Organize the experiments: summary across recordings
load(fullfile(save_to,'recList.mat'))
n=height(recList);
recList.Ncells=nan(n,1); recList.Nactive=nan(n,1); recList.MedianRate=nan(n,1);
recList.NCS=nan(n,1); recList.CSrate=nan(n,1); recList.ResultFile=repmat({''},n,1);
if ~ismember('Extracted',recList.Properties.VariableNames), recList.Extracted=repmat({''},n,1); end
recList.Analyzed=repmat({''},n,1);
for f=1:n
    resFile=out_file(recList.Path{f},'Result.mat');
    if ~exist(resFile,'file'), continue; end
    try, R=load(resFile); R=R.Result; catch, fprintf(2,'  UNREADABLE %s\n',resFile); continue; end
    recList.ResultFile{f}=resFile;
    if isfield(R,'extractedOn'), recList.Extracted{f}=R.extractedOn; end   % else keep the recorded one
    if isfield(R,'analyzedOn'), recList.Analyzed{f}=R.analyzedOn;
    else, d=dir(resFile); recList.Analyzed{f}=char(datetime(d.datenum,'ConvertFrom','datenum','Format','yyyy-MM-dd HH:mm'));
    end
    recList.Ncells(f)=size(R.ftprnt,3);
    if isfield(R,'isActive')
        recList.Nactive(f)=sum(R.isActive);
        recList.MedianRate(f)=median(R.spRate(R.isActive));
        recList.NCS(f)=sum(R.nCS);
        recList.CSrate(f)=median(R.nCS(R.isActive))/(sum(~R.excluded)/R.frm_rate);
    end
end
write_rec_sheet(recList,xlsFile,recSheet);
save(fullfile(save_to,'recList.mat'),'recList');

figure(400); clf; theme(gcf,'light');
set(gcf,'Color','w','Position',[100 100 1400 800],'Visible',figVis);
tiledlayout(2,2,'TileSpacing','compact');
prots=unique(recList.Protocol,'stable');
cnt=zeros(numel(prots),2);
for p=1:numel(prots)
    cnt(p,1)=sum(strcmp(recList.Protocol,prots{p}) & recList.Use==1);
    cnt(p,2)=sum(strcmp(recList.Protocol,prots{p}) & recList.Use==0);
end
nexttile; barh(categorical(prots,prots),cnt,'stacked');
legend({'motion corrected','no mc.bin (excluded)'},'Location','southeast');
xlabel('# recordings'); title('Protocols used');
done=~isnan(recList.Ncells);
nexttile; hold on
for p=1:numel(prots)
    ii=find(done & strcmp(recList.Protocol,prots{p}));
    scatter(p+0.1*randn(numel(ii),1),recList.Nactive(ii),25,'filled');
end
set(gca,'xtick',1:numel(prots),'xticklabel',prots); xtickangle(30);
ylabel('# active cells / recording'); title('Active cells');
nexttile; hold on
for p=1:numel(prots)
    ii=find(done & strcmp(recList.Protocol,prots{p}));
    scatter(p+0.1*randn(numel(ii),1),recList.MedianRate(ii),25,'filled');
end
set(gca,'xtick',1:numel(prots),'xticklabel',prots,'yscale','log'); xtickangle(30);
ylabel('median firing rate (Hz)'); title('Firing rate');
nexttile; hold on
for p=1:numel(prots)
    ii=find(done & strcmp(recList.Protocol,prots{p}));
    scatter(p+0.1*randn(numel(ii),1),recList.CSrate(ii),25,'filled');
end
set(gca,'xtick',1:numel(prots),'xticklabel',prots); xtickangle(30);
ylabel('median CS rate, active cells (Hz)'); title('Complex spikes');
exportgraphics(gcf,fullfile(save_to,'Summary_protocols.png'),'Resolution',150);

%% 5. Representative movie: reconstructed from footprints x traces
% each neuron in its own color over the gray mean image, brighter with depolarization
movRec='BHLpop081';    % Rec ID
movWin=[3.3 9.6];      % s
M=struct('bin',2,'fps',60,'smooth',2,'zLow',1,'zHigh',8,'baseAlpha',0.25, ...
    'bgGain',0.5,'activeOnly',true);
load(fullfile(save_to,'recList.mat'))
f=find(strcmp(recList.RecID,movRec));
load(out_file(recList.Path{f},'Result.mat'))
frm=max(1,round(movWin(1)*Result.frm_rate)):min(Result.nFrames,round(movWin(2)*Result.frm_rate));
make_recon_movie(Result,frm,out_file(Result.fpath,'recon.mp4'),M);

%% Functions
function fn=out_file(fpath,suffix)
% per-recording outputs live next to the voltage movie
fn=fullfile(fpath,['PopReanalysis_' suffix]);
end

function recDirs=find_recordings(parent)
L=[dir(fullfile(parent,'**','settings.mat')); dir(fullfile(parent,'**','Sq_camera.bin')); ...
    dir(fullfile(parent,'**','mc*.bin'))];
L=L(~startsWith({L.name},'._'));
recDirs=unique({L.folder})';
skip=contains(recDirs,{'hilo','snaps'},'IgnoreCase',true) | ...
    ~cellfun(@isempty,regexp(recDirs,'_y\d+_z\d+$','once'));
recDirs=recDirs(~skip);
end

function rec=inspect_recording(d,parent)
[~,nm]=fileparts(d);
[~,pnm]=fileparts(parent);
rec.RecID='';
rec.Date=str2double(regexp(pnm,'^\d{8}','match','once'));
rec.Path=d;
rec.Parent=parent;
rec.Recording=nm;
rec.Relative=strrep(d,[parent '/'],'');

tok=regexp(nm,'BHLm\d+','match','once');
if isempty(tok), tok=regexp(nm,'mouse\d+','match','once'); end
if isempty(tok), tok=regexp(pnm,'BHLm\d+','match','once'); end
rec.Mouse=tok;
tok=regexp(lower(rec.Relative),'(?:^|[_/])(\d+x)(?=[_/]|$)','tokens');
if isempty(tok), rec.Mag=''; else, rec.Mag=tok{end}{1}; end
tok=regexp(nm,'FOV\d+','match','once'); rec.FOV=tok;
rec.BlueName=parse_blue(nm);

% rig settings
sf=dir(fullfile(d,'*.set')); sf=sf(~startsWith({sf.name},'._'));
rec.Waveform=''; camFr=NaN; expo=NaN;
if ~isempty(sf)
    txt=fileread(fullfile(d,sf(1).name));
    rec.Waveform=get_set_field(txt,'WaveformFunctionDropDown');
    camFr=str2double(get_set_field(txt,'camFramesInEditField'));
    expo=str2double(get_set_field(txt,'ExposureTimeMsEditField'));
end
[W,H]=read_frame_size(d);
rec.NDMDpatt=0;
sfile=fullfile(d,'settings.mat');
if exist(sfile,'file')
    try
        w=whos('-file',sfile);
        if ismember('dmd_mask_sequence_rois',{w.name})
            S=load(sfile,'dmd_mask_sequence_rois');
            if iscell(S.dmd_mask_sequence_rois), rec.NDMDpatt=numel(S.dmd_mask_sequence_rois); end
        end
    catch
        rec.NDMDpatt=NaN;   % unreadable settings.mat
    end
end
csvf=dir(fullfile(d,'*.csv')); csvf=csvf(~startsWith({csvf.name},'._'));
rec.TreadmillCsv=~isempty(csvf);

% motion correction status
mcf=dir(fullfile(d,'mc*.bin')); mcf=mcf(~startsWith({mcf.name},'._'));
isChunk=~cellfun(@isempty,regexp({mcf.name},'^mc\d+\.bin$','once'));
hasSingle=any(strcmp({mcf.name},'mc.bin'));
if hasSingle
    rec.MCstatus='mc.bin';
    hasTrace=exist(fullfile(d,'mcTrace.mat'),'file')>0;
elseif any(isChunk)
    rec.MCstatus=sprintf('mc01-%02d.bin',sum(isChunk));
    hasTrace=exist(fullfile(d,'mcTrace01.mat'),'file')>0;
else
    rec.MCstatus='NO MC'; hasTrace=false;
end
rec.NFrames=camFr;
if (hasSingle || any(isChunk)) && ~isnan(W)
    ch=get_mc_chunks(d,W,H); rec.NFrames=ch(end).gStart+ch(end).nUse-1;
elseif exist(fullfile(d,'Sq_camera.bin'),'file') && ~isnan(W)
    s=dir(fullfile(d,'Sq_camera.bin')); rec.NFrames=s.bytes/(W*H*2);
end
rec.FrameRate=1000/(expo+0.02);
rec.Duration=rec.NFrames/rec.FrameRate;

[rec.Protocol,rec.ProtocolDetail,rec.State]=classify_protocol(nm,rec.Waveform, ...
    rec.NDMDpatt,rec.TreadmillCsv,rec.NFrames,expo);
rec.Use=double(~strcmp(rec.MCstatus,'NO MC') && hasTrace && ~isnan(W));
if strcmp(rec.MCstatus,'NO MC')
    rec.Misc='No motion-corrected movie (mc.bin): likely bad recording, excluded';
elseif ~hasTrace
    rec.Misc='mc.bin without mcTrace.mat, excluded';
else
    rec.Misc='';
end
end

function [prot,detail,state]=classify_protocol(nm,wf,nPatt,hasCsv,nFr,expo)
n=lower(nm);
tk=strsplit(n,{'_','-'});
state='';
if contains(n,'awake'), state='Awake';
elseif any(ismember(tk,{'ane','anes'})), state='Anesthetized'; end
freq=regexp(n,'\d+hz','match','once');
if expo>=20
    prot='Slow imaging'; detail=sprintf('%g ms exposure (not voltage rate)',expo);
elseif contains(n,'treadmill') && (hasCsv || nFr>50000)
    prot='Treadmill'; state='Awake';
    if contains(n,'const'), detail='Treadmill + constant blue';
    else, detail='Treadmill, spontaneous'; end
elseif contains(n,'compete')
    prot='Selective stim';
    detail=sprintf('Assembly competition, %d DMD patterns (%s)',nPatt,regexp(n,'compete\d+v\d+','match','once'));
elseif contains(n,'cluster')
    prot='Selective stim';
    detail=sprintf('Cluster stimulation, %d DMD patterns (%s)',nPatt,regexp(n,'cluster\d','match','once'));
elseif any(ismember(tk,{'bis'})) || contains(n,'bistab')
    prot='Bistability'; detail=['2-color, waveform ' wf];
elseif contains(n,{'sqr_train','pulse_train'})
    prot='Pulse train'; detail=['Blue pulse train ' freq];
elseif contains(n,{'const','contst'})
    prot='Constant stim'; detail='Constant blue';
elseif any(ismember(tk,{'rp','pr','bp'})) || contains(n,'pulse_ramp')
    prot='Pulse + ramp'; detail='Wide-field blue pulses + ramp';
elseif contains(n,{'pulse','blue_p','bluep'})
    prot='Pulse stim'; detail='Wide-field blue pulses';
elseif contains(n,'treadmill')
    prot='Treadmill test'; detail='Short treadmill test, no position data';
elseif contains(wf,'SimpleWaveform') && nFr==8000
    prot='Pulse + ramp'; detail='Inferred from waveform (8000 frames, simple waveform)';
else
    prot='Other'; detail='';
end
end

function b=parse_blue(nm)
b=NaN;
tok=regexp(nm,'(?<=_)B(\d*)p(\d+)','tokens','once');
if ~isempty(tok)
    if isempty(tok{1}), tok{1}='0'; end
    b=str2double([tok{1} '.' tok{2}]); return
end
tok=regexp(nm,'(?<=_)B(\d+)(?=_|$)','tokens','once');
if ~isempty(tok), b=str2double(tok{1}); return; end
tok=regexp(nm,'[Bb]lue_?[Pp](\d+)','tokens','once');
if ~isempty(tok), b=str2double(['0.' tok{1}]); end
end

function [W,H]=read_frame_size(d)
W=NaN; H=NaN;
ep=fullfile(d,'experimental_parameters.txt');
if ~exist(ep,'file'), return; end
tok=regexp(fileread(ep),'pixel:_*\s*(\d+)','tokens');
if numel(tok)>=2, W=str2double(tok{1}{1}); H=str2double(tok{2}{1}); end
end

function v=get_set_field(txt,key)
tok=regexp(txt,[key ':\s*"?([^"\r\n]*)"?'],'tokens','once');
if isempty(tok), v=''; else, v=strtrim(tok{1}); end
end

function write_rec_sheet(recList,xlsFile,recSheet)
hdr={'Rec ID','Date','Data Location','Parent folder','Recording','Mouse#','Mag','FOV', ...
    'Protocol','Protocol detail','Awake/Anesthetized','Blue (from name)','Waveform', ...
    'N DMD patt','N frames','Frame rate (Hz)','Duration (s)','Treadmill csv','MC status','Use?','Misc'};
C=[recList.RecID num2cell(recList.Date) recList.Path recList.Parent recList.Recording ...
    recList.Mouse recList.Mag recList.FOV recList.Protocol recList.ProtocolDetail recList.State ...
    num2cell(recList.BlueName) recList.Waveform num2cell(recList.NDMDpatt) ...
    num2cell(recList.NFrames) num2cell(round(recList.FrameRate,1)) num2cell(round(recList.Duration,1)) ...
    num2cell(double(recList.TreadmillCsv)) recList.MCstatus num2cell(recList.Use) recList.Misc];
if ismember('Ncells',recList.Properties.VariableNames)
    hdr=[hdr {'N cells','N active','Median rate (Hz)','N complex spikes','Median CS rate (Hz)', ...
        'Result file','Extracted','Analyzed'}];
    C=[C num2cell(recList.Ncells) num2cell(recList.Nactive) num2cell(round(recList.MedianRate,2)) ...
        num2cell(recList.NCS) num2cell(round(recList.CSrate,3)) recList.ResultFile recList.Extracted recList.Analyzed];
end
C(cellfun(@(x) isnumeric(x) && isscalar(x) && isnan(x),C))={[]};
writecell([{'Population data, one row per recording. Use?=0: excluded (see Misc).'}; ...
    {['Updated ' char(datetime('now','Format','yyyy-MM-dd HH:mm'))]}],xlsFile, ...
    'Sheet',recSheet,'WriteMode','overwritesheet');
writecell([hdr; C],xlsFile,'Sheet',recSheet,'Range','A4');
end

function ch=get_mc_chunks(d,W,H)
% mc.bin, or mc01..mcNN.bin written in 10000-frame steps with 11 overlapping frames
fb=W*H*2;
if exist(fullfile(d,'mc.bin'),'file')
    s=dir(fullfile(d,'mc.bin'));
    ch=struct('file',fullfile(d,'mc.bin'),'mcFile',fullfile(d,'mcTrace.mat'), ...
        'nUse',floor(s.bytes/fb),'gStart',1);
    return
end
L=dir(fullfile(d,'mc*.bin'));
L=L(~cellfun(@isempty,regexp({L.name},'^mc\d+\.bin$','once')));
[~,o]=sort({L.name}); L=L(o);
nFile=floor([L.bytes]/fb);
step=floor(nFile(1)/1000)*1000;
ch=struct('file',{},'mcFile',{},'nUse',{},'gStart',{});
for k=1:numel(L)
    ch(k).file=fullfile(d,L(k).name);
    ch(k).mcFile=fullfile(d,sprintf('mcTrace%02d.mat',k));
    ch(k).nUse=min(nFile(k),step);
    if k==numel(L), ch(k).nUse=nFile(k); end
    ch(k).gStart=(k-1)*step+1;
end
end

function [mov,mcT]=read_chunk(ch,W,H)
fid=fopen(ch.file,'r');
mov=fread(fid,W*H*ch.nUse,'*uint16',0,'l');
fclose(fid);
mov=permute(reshape(mov,W,H,ch.nUse),[2 1 3]);
M=load(ch.mcFile,'mcTrace');
mcT=double(M.mcTrace(1:min(end,ch.nUse),:));
if size(mcT,1)<ch.nUse, mcT(end+1:ch.nUse,:)=repmat(mcT(end,:),ch.nUse-size(mcT,1),1); end
end

function C=correct_chunk(mov,mcT,blue,P)
% bleaching: divide by an exponential fit of the FOV mean (blue-off frames)
% motion: regress x, y, x^2, y^2, xy shifts and x, y velocity out of every pixel, block-wise
% frames with the imaging laser off or fast motion transients are replaced by the
% block mean image (they are excluded later)
mov=single(mov); T=size(mov,3);
mF=double(squeeze(mean(mov,[1 2])));
invalid=mF<0.5*median(mF);
invalid=imdilate(invalid,ones(5,1));
blue=double(blue(:));
fitFr=~invalid;
if any(blue>0)
    [~,bo]=get_blueoffTrace(mF',blue',70);
    if sum(bo(:) & fitFr)>0.1*T, fitFr=bo(:) & fitFr; end
end
x=find(fitFr);
try
    bleach=expfitDM_2(x,mF(x),(1:T)',1000);
    if any(~isfinite(bleach) | bleach<=0), error('bad fit'); end
catch
    bleach=interp1(x,movmean(mF(x),500),(1:T)','linear','extrap');
end
mov=mov./reshape(single(bleach),1,1,[]);

% fast motion transients: robust z of the frame-to-frame RMS change, blue on/off edges spared
dd=zeros(T,1);
for b=1:1000:T
    bj=max(b-1,1):min(T,b+999);
    dd(bj(2:end))=squeeze(sqrt(mean(diff(mov(:,:,bj),1,3).^2,[1 2])));
end
dd(1)=median(dd);
dz=(dd-median(dd))/(1.4826*mad(dd,1)+eps);
artifact=dz>P.artZ;
artifact(imdilate([false; diff(blue>0)~=0],ones(7,1)))=false;
artifact=imdilate(artifact,ones(2*P.artDilate+1,1)) & ~invalid;
bad=invalid | artifact;

vel=[zeros(1,2); diff(mcT,1,1)];
nB=ceil(T/P.blockLen);
meanImg=zeros(size(mov,1),size(mov,2),nB,'single');
for b=1:nB
    bi=(b-1)*P.blockLen+1:min(T,b*P.blockLen);
    blk=mov(:,:,bi);
    ok=~bad(bi);
    if ~any(ok), ok(:)=true; end
    meanImg(:,:,b)=mean(blk(:,:,ok),3);
    if any(~ok), blk(:,:,~ok)=repmat(meanImg(:,:,b),1,1,sum(~ok)); end
    m=mcT(bi,:); m=m-mean(m,1); v=vel(bi,:);
    mov(:,:,bi)=SeeResiduals(blk,[m m.^2 m(:,1).*m(:,2) v]');
end
C.res=mov; C.meanImg=meanImg; C.bleach=bleach; C.invalid=invalid; C.artifact=artifact;
C.mcT=mcT; C.mF=mF;
end

function Result=reextract_recording(rec,P,interactiveUI)
rng(0,'twister');   % kmeans in extractVoltronST: the same footprints on every run
fpath=rec.Path{1};
[sz1,sz2]=read_frame_size(fpath);   % horizontal, vertical
S=load(fullfile(fpath,'settings.mat'),'DAQ_waves','dmd_mask_sequence_rois');
chunks=get_mc_chunks(fpath,sz1,sz2);
nFrames=chunks(end).gStart+chunks(end).nUse-1;

amp=S.DAQ_waves.amplitude;
interval=median(diff(find(diff(amp(1,1:min(5000,size(amp,2))))>0)));   % DAQ samples (10 us) per frame
frmT=interval*1e-5;
idx=min(round((1:nFrames)*interval),size(amp,2));

Result=struct;
Result.fpath=fpath; Result.RecID=rec.RecID{1}; Result.recName=rec.Recording{1};
Result.protocol=rec.Protocol{1}; Result.protocolDetail=rec.ProtocolDetail{1};
Result.state=rec.State{1};
Result.frm_rate=1/frmT; Result.nFrames=nFrames; Result.P=P;
Result.DMDPatt={}; Result.clist_p={};
if iscell(S.dmd_mask_sequence_rois) && numel(S.dmd_mask_sequence_rois)>1
    Result.DMDPatt=S.dmd_mask_sequence_rois;
end
[Result.Blue,Result.Orange,Result.BluePatt]=get_stim(amp,interval,nFrames,numel(Result.DMDPatt));
Result.stimVer=2;
clear S amp idx

Result.Arduino=[];
csvf=dir(fullfile(fpath,'*.csv')); csvf=csvf(~startsWith({csvf.name},'._'));
if ~isempty(csvf)
    try
        D=load(fullfile(fpath,'settings.mat'),'DAQ_data');
        A=match_treadmill_DAQ(fullfile(fpath,csvf(1).name),1e-5,D.DAQ_data,frmT,2); % time, position, reward, run
        A=fitlen(A',nFrames)'; A(:,1)=(1:nFrames)'*frmT;
        Result.Arduino=A;
    catch ME
        warning('Treadmill matching failed: %s',ME.message);
    end
    clear D
end

% detection windows, each inside one chunk
nDet=min(ceil(nFrames/P.detWinLen),P.maxDetWin);
if nFrames>20000, nDet=min(ceil(nFrames/P.detWinLen),P.maxDetWinLong); end
edges=round(linspace(0,nFrames,nDet+1));
detWin=zeros(nDet,2);
for k=1:nDet
    L=min(P.detWinLen,edges(k+1)-edges(k));
    s=round((edges(k)+edges(k+1))/2-L/2)+1;
    c=find([chunks.gStart]<=max(s,1),1,'last');
    L=min(L,chunks(c).nUse);
    cEnd=chunks(c).gStart+chunks(c).nUse-1;
    s=max(min(s,cEnd-L+1),chunks(c).gStart);
    detWin(k,:)=[s s+L-1];
end

vo=P.vopts; vo.doPlot=false; vo.uiSeedThresh=false; vo.uiClusterSize=false; vo.uiSynFilter=false;
vo.synFilter=struct('ampRange',[-inf inf],'countRange',[vo.minEventsSyn inf], ...
    'areaRange',[vo.minArea P.maxArea]);
edgeMsk=false(sz2,sz1);
edgeMsk(P.edgeMask+1:end-P.edgeMask,P.edgeMask+1:end-P.edgeMask)=true;
tileRes={}; cPrev=0; Result.detSkipped=false(1,nDet);
for k=1:nDet
    c=find([chunks.gStart]<=detWin(k,1),1,'last');
    if c~=cPrev
        [mov,mcT]=read_chunk(chunks(c),sz1,sz2);
        g=chunks(c).gStart+(0:chunks(c).nUse-1);
        CC=correct_chunk(mov,mcT,Result.Blue(g),P); clear mov
        cPrev=c;
    end
    wi=(detWin(k,1):detWin(k,2))-chunks(c).gStart+1;
    vo_k=vo;
    if interactiveUI && k==1
        vo_k.doPlot=true; vo_k.uiSeedThresh=true; vo_k.uiClusterSize=true; vo_k.uiSynFilter=true;
        vo_k.diagSaveDir=out_file(fpath,'diag'); vo_k.diagPrefix=Result.RecID;
    end
    try
        R=extractVoltronST(movmean(-CC.res(:,:,wi),P.detSmooth,3).*edgeMsk,mean(CC.meanImg,3),vo_k);
    catch ME   % e.g. an artifact episode with >1e5 events (pdist2 out of memory)
        fprintf(2,'  detection window %d skipped: %s\n',k,ME.message);
        Result.detSkipped(k)=true;
        continue
    end
    if interactiveUI && k==1   % carry the chosen thresholds to the other windows
        vo.zSeedMin=R.opts.zSeedMin; vo.r_max=R.opts.r_max; vo.synFilter=R.opts.synFilter;
    end
    if ~isempty(R.S_glu)
        tileRes{end+1}=struct('S_glu',R.S_glu,'synCentroid',R.synCentroid);
    end
    fprintf('  detection window %d/%d [%d-%d]: %d footprints\n',k,nDet, ...
        detWin(k,1),detWin(k,2),size(R.S_glu,3)*~isempty(R.S_glu));
end
Result.detWin=detWin;
if isempty(tileRes)
    Result.ftprnt=[]; Result.centroid=[]; Result.ref_im=mean(CC.meanImg,3);
    return
end
Dict=mergeSynapseDict(tileRes,P.mergeRadius);
Sglu=double(Dict.S_glu);
nCell=size(Sglu,3);
Result.ftprnt=single(Sglu); Result.centroid=Dict.centroid;

% projection over the whole recording (= extractVoltronST project mode, footprint pixels only)
Smat=reshape(Sglu,[],nCell);
usedPix=find(max(Smat,[],2)>0);
pS=single(pinv(Smat(usedPix,:)));
Result.traces=zeros(nCell,nFrames,'single'); Result.F0=zeros(nCell,nFrames,'single');
Result.mcTrace=zeros(nFrames,2); Result.bleach=zeros(nFrames,1);
Result.meanF=zeros(nFrames,1); Result.invalid=false(1,nFrames); Result.artifact=false(1,nFrames);
for c=1:numel(chunks)
    g=chunks(c).gStart+(0:chunks(c).nUse-1);
    if c~=cPrev
        [mov,mcT]=read_chunk(chunks(c),sz1,sz2);
        CC=correct_chunk(mov,mcT,Result.Blue(g),P); clear mov
        cPrev=c;
    end
    if c==1, Result.ref_im=mean(CC.meanImg,3); end
    M=reshape(CC.res,[],numel(g));
    for b=1:ceil(numel(g)/P.blockLen)
        bi=(b-1)*P.blockLen+1:min(numel(g),b*P.blockLen);
        Result.traces(:,g(bi))=-pS*M(usedPix,bi);
        img=CC.meanImg(:,:,b);
        Result.F0(:,g(bi))=repmat(pS*img(usedPix),1,numel(bi));
    end
    clear M
    Result.mcTrace(g,:)=CC.mcT; Result.bleach(g)=CC.bleach;
    Result.meanF(g)=CC.mF; Result.invalid(g)=CC.invalid; Result.artifact(g)=CC.artifact;
    fprintf('  chunk %d/%d projected\n',c,numel(chunks));
end
Result.dFF=Result.traces./Result.F0;

for j=1:numel(Result.DMDPatt)   % cells inside each DMD pattern
    roi=dmd_pat2roi(Result.DMDPatt{j},[sz2 sz1]);
    cc=round(Result.centroid);
    cc(:,1)=min(max(cc(:,1),1),sz1); cc(:,2)=min(max(cc(:,2),1),sz2);
    Result.clist_p{j}=find(roi(sub2ind([sz2 sz1],cc(:,2),cc(:,1))));
end
Result.extractedOn=char(datetime('now','Format','yyyy-MM-dd HH:mm'));
end

function [Blue,Orange,BluePatt]=get_stim(amp,interval,nFrames,nP)
% blue: max of the DAQ waveform within each frame (pulses shorter than a frame are kept)
% BluePatt: same per DMD pattern (pattern index from the DMD trigger, as show_multiROI_waveform)
nF=min(nFrames,floor(size(amp,2)/interval));
k=reshape(1:nF*interval,interval,nF);
b=reshape(amp(4,k),interval,nF);
Blue=fitlen(max(b,[],1),nFrames);
Orange=fitlen(mean(reshape(amp(5,k),interval,nF),1),nFrames);
BluePatt=[];
if nP>1
    dmd_time=mod(cumsum(amp(3,:)),nP); dmd_time(dmd_time==0)=nP;
    dt=reshape(dmd_time(k),interval,nF);
    BluePatt=zeros(nP,nFrames);
    for j=1:nP
        BluePatt(j,:)=fitlen(max(b.*(dt==j),[],1),nFrames);
    end
end
end

function Result=refresh_stim(Result)
% results extracted before the per-frame max sampling: recompute the stimulus
if isfield(Result,'stimVer') && Result.stimVer>=2, return; end
S=load(fullfile(Result.fpath,'settings.mat'),'DAQ_waves');
amp=S.DAQ_waves.amplitude;
interval=median(diff(find(diff(amp(1,1:min(5000,size(amp,2))))>0)));
[Result.Blue,Result.Orange,Result.BluePatt]=get_stim(amp,interval,Result.nFrames,numel(Result.DMDPatt));
Result.stimVer=2;
end

function safe_save(resFile,Result)
% write to a temporary file, check that it loads back completely, then replace the
% result file (a bad write never overwrites a good result)
tmp=[resFile(1:end-4) '_tmp.mat'];
for attempt=1:2
    save(tmp,'Result','-v7.3');
    try
        S=load(tmp); clear S
        movefile(tmp,resFile,'f');
        return
    catch ME
        fprintf(2,'  save check failed (%s), attempt %d\n',ME.message,attempt);
        if exist(tmp,'file'), delete(tmp); end
    end
end
error('Could not write a readable result file: %s',resFile);
end

function x=fitlen(x,n)
if size(x,2)>=n, x=x(:,1:n); else, x(:,end+1:n)=repmat(x(:,end),1,n-size(x,2)); end
end

function Result=analyze_recording(Result,P,csUI)
rng(0,'twister');   % surrogates and null samples: the same result on every run
Result.analyzedOn=char(datetime('now','Format','yyyy-MM-dd HH:mm'));
tr=double(Result.dFF); tr(~isfinite(tr))=0;
valid=~(Result.invalid | Result.artifact);
hi=tr-movmedian(tr,P.hiWin,2);
noise=get_threshold(hi,1);
Result.zTr=single(hi./noise);
% motion-like frames: many cells deflect strongly negative together (residual motion swings
% both ways, depolarizing population events do not), so synchronous events are kept
Result.globalArt=false(1,Result.nFrames);
if size(tr,1)>=10
    fn=double(mean(Result.zTr<-P.negZ,1));
    thr=max(P.globalFrac,median(fn)+6*1.4826*mad(fn,1));
    Result.globalArt=imdilate(fn>thr,ones(1,2*P.globalDilate+1));
end
sp=find_spike_bh(double(Result.zTr),P.spThr,P.spProm);
valid=valid & ~Result.globalArt;
Result.excluded=~valid;
sp(:,~valid)=0;
Result.spike=sp>0;
Result=detect_cs(Result,P,csUI);   % also adds / removes spikes on complex-spike plateaus
Result.nSpikes=sum(Result.spike,2);
Result.spRate=Result.nSpikes/(sum(valid)/Result.frm_rate);
Result.isActive=Result.nSpikes>=P.minSpikes;
Result.placeMap=[];
if ~isempty(Result.Arduino)
    Result.placeMap=place_tuning(Result.spike,Result.Arduino,valid,Result.frm_rate,P.placeBin);
end
end

function make_recon_movie(Result,frm,outFile,M)
% reconstructed movie = sum_n footprint_n x color_n x activity_n(t) on the gray mean image
% activity: z-scored trace, 0 at M.zLow and 1 at M.zHigh, each neuron at M.baseAlpha at rest
cells=1:size(Result.ftprnt,3);
if M.activeOnly && isfield(Result,'isActive'), cells=find(Result.isActive)'; end
S=double(Result.ftprnt(:,:,cells));
S=S./max(S,[],[1 2]);
[H,W,nC]=size(S);
Sv=reshape(S,H*W,nC);
rng(1); cmap=hsv(nC); cmap=cmap(randperm(nC),:);   % neighbours get different hues

z=movmean(double(Result.zTr(cells,:)),M.smooth,2);
a=min(max((z-M.zLow)/(M.zHigh-M.zLow),0),1);
a(:,Result.excluded)=0;
a=M.baseAlpha+(1-M.baseAlpha)*a;
nb=floor(numel(frm)/M.bin);
fb=reshape(frm(1:nb*M.bin),M.bin,nb);
aB=squeeze(max(reshape(a(:,fb(:)),nC,M.bin,nb),[],2));   % max within each bin keeps spikes

ref=double(Result.ref_im);
bg=mat2gray(ref,double(prctile(ref(:),[1 99.7])))*M.bgGain;
bg=bg(:);
t=(1:Result.nFrames)/Result.frm_rate;

fig=figure('Visible','off','Color','k','Position',[50 50 1000 760]); theme(fig,'dark');
ax1=axes(fig,'Position',[0.02 0.25 0.96 0.72]);
hIm=image(ax1,zeros(H,W,3)); axis(ax1,'image','off');
hT=text(ax1,5,12,'','Color','w','FontSize',14,'FontWeight','bold');
title(ax1,sprintf('%s  %s  (%d neurons)',Result.RecID,Result.recName,nC), ...
    'Color','w','Interpreter','none','FontSize',11);
ax2=axes(fig,'Position',[0.06 0.05 0.9 0.15]);
if ~isempty(Result.BluePatt), stim=Result.BluePatt(any(Result.BluePatt,2),:); else, stim=Result.Blue; end
plot(ax2,t(frm),stim(:,frm)','LineWidth',1); hold(ax2,'on');
hC=xline(ax2,t(frm(1)),'r','LineWidth',1.5);
axis(ax2,'tight'); xlabel(ax2,'time (s)'); ylabel(ax2,'blue (V)');

v=VideoWriter(outFile,'MPEG-4'); v.FrameRate=M.fps; v.Quality=95;
open(v);
for k=1:nb
    al=Sv.*aB(:,k)';
    alpha=min(sum(al,2),1);
    rgb=bg.*(1-alpha)+min(al*cmap,1);
    hIm.CData=reshape(min(rgb,1),H,W,3);
    hT.String=sprintf('%.3f s',t(fb(1,k)));
    hC.Value=t(fb(1,k));
    drawnow;
    writeVideo(v,getframe(fig));
end
close(v); close(fig);
fprintf('Movie saved: %s (%d frames, %.1fx slower than real time)\n',outFile,nb, ...
    Result.frm_rate/(M.fps*M.bin));
end

function [z,sub]=cs_signals(Result,P)
% z: dF/F in units of the spike-detection noise (the trace given to interactive_cs_finder)
% sub: its envelope as computed inside interactive_cs_finder (movmean 20 fr - movmedian 500 fr)
tr=double(Result.dFF); tr(~isfinite(tr))=0;
z=(tr-median(tr,2))./get_threshold(tr-movmedian(tr,P.hiWin,2),1);
sub=movmean(z,20,2); sub=sub-movmedian(sub,500,2);
end

function o=cs_opts(P)
o=struct('CS_thres',P.cs.thres,'N_Spike',P.cs.nSpike,'Max_AUC',P.cs.maxAUC, ...
    'N_Spike2ISI',P.cs.nSpikeISI,'showUI',false);
end

function Result=detect_cs(Result,P,csUI)
% complex spikes with interactive_cs_finder; with csUI its window opens on the most active
% cell and the settings chosen there are applied to every cell.
% Before that, spikes inside putative CS are re-detected on the plateau (cs_spikes); the CS
% onset is the first spike of each CS, always a regularly detected one.
[z,sub]=cs_signals(Result,P);
sp=double(Result.spike);
o=cs_opts(P);
[spAll,added,removed]=cs_spikes(z,sub,sp,~Result.excluded,o,P);
if csUI && any(sum(sp,2)>=P.minSpikes)
    [~,n]=max(sum(sp,2));
    o.showUI=true;
    [~,o.CS_thres,prm]=interactive_cs_finder(spAll(n,:),z(n,:),o);
    o.N_Spike=prm.N_Spike; o.Max_AUC=prm.Max_AUC; o.N_Spike2ISI=prm.N_Spike2ISI;
    o.showUI=false;
    [spAll,added,removed]=cs_spikes(z,sub,sp,~Result.excluded,o,P);   % regions follow the new thresholds
end
CS=false(size(z)); on=false(size(z));
for n=find(sum(spAll,2)>=o.N_Spike)'
    cs=interactive_cs_finder(spAll(n,:),z(n,:),o)>0;
    lab=bwlabel(cs);
    for b=1:max(lab)
        idx=find(lab==b);
        f1=idx(find(sp(n,idx)>0 & ~removed(n,idx),1));   % CS onset = first regularly detected spike
        if isempty(f1) || sum(spAll(n,f1:idx(end)))<o.N_Spike, cs(idx)=false; continue; end
        cs(idx(1):f1-5)=false;         % span starts 4 frames before it, as in the finder
        on(n,f1)=true;
    end
    CS(n,:)=cs;
end
added=added & CS; removed=removed & CS;   % re-detection is kept only inside a detected CS
Result.spike=(Result.spike & ~removed) | added;
Result.spikeAdded=added; Result.spikeRemoved=removed;
Result.CStrace=CS; Result.CSonset=on; Result.nCS=sum(on,2);
Result.CSparams=rmfield(o,'showUI');
% one row per CS: cell, first-spike frame, first-spike time (s, sub-frame peak), last frame,
% number of spikes
[c,fr]=find(on); [c,o2]=sort(c); fr=fr(o2);
Result.CSlist=zeros(numel(c),5);
for k=1:numel(c)
    y=z(c(k),max(1,fr(k)-1):min(end,fr(k)+1)); d=0;
    if numel(y)==3 && y(1)-2*y(2)+y(3)<0, d=0.5*(y(1)-y(3))/(y(1)-2*y(2)+y(3)); end
    e=fr(k)-2+find(~CS(c(k),fr(k):end),1); if isempty(e), e=Result.nFrames; end
    Result.CSlist(k,:)=[c(k) fr(k) (fr(k)+d)/Result.frm_rate e sum(Result.spike(c(k),fr(k):e))];
end
end

function [spAll,added,removed]=cs_spikes(z,sub,sp,valid,o,P)
% spikes inside putative CS = envelope transient (above CS_thres(2), reaching CS_thres(1)) plus
% the 20 frames before it (spikes the finder joins), after the first regular spike, on z
% without removing the plateau. Added: peaks > 2 frames from any regular spike, reaching
% P.spThr above the baseline in front of the plateau, with a 7-pt amplitude >= P.cs.prom.
% On plateaus >= noisyRatio x noisier than quiet baseline nothing is added and regular spikes
% that do not stand out from that noise are removed.
added=false(size(sp)); removed=false(size(sp));
for n=find(any(sp,2))'
    L=bwlabel(sub(n,:)>o.CS_thres(2));
    k=unique(L(movmean(sub(n,:),3)>o.CS_thres(1))); k(k==0)=[];
    if isempty(k), continue; end
    R=bwlabel(movmax(double(ismember(L,k)),[0 20])>0);
    base=nan(size(R));   % baseline in front of the plateau, on the frames after the first spike
    for r=1:max(R)
        ii=find(R==r); f1=ii(find(sp(n,ii)>0,1));
        if isempty(f1) || f1<4, continue; end
        base(f1+1:ii(end))=median(z(n,max(1,f1-2-P.cs.baseFr):f1-3));
    end
    quiet=sub(n,:)<o.CS_thres(2) & ~(movmax(sp(n,:),[20 20])>0) & ...
        ~(movmax(double(R>0),[30 30])>0) & valid;
    nearSp=movmax(sp(n,:),[2 2])>0;
    d=[diff(z(n,:)) 0]; q=quiet & [quiet(2:end) false];
    if ~any(q), q=true(size(d)); end
    sq=1.4826*mad(d(q),1)/sqrt(2);   % quiet-baseline noise, from frame-to-frame differences
    d(nearSp | [nearSp(2:end) false])=NaN;
    ratio=1.4826*movmad(d,P.cs.noiseWinFr,'omitnan')/sqrt(2)/sq;   % local noise / quiet noise
    ratio(~isfinite(ratio))=1;
    noisy=ratio>=P.cs.noisyRatio;
    a7=z(n,:)-movmedian(z(n,:),7);   % 7-pt amplitude, the prominence of find_spike_bh
    on=~isnan(base) & valid;
    removed(n,:)=sp(n,:)>0 & on & noisy & a7<P.spProm*ratio;
    [~,fr]=findpeaks(z(n,:),'MinPeakDistance',P.cs.minIsiFr);
    fr=fr(on(fr) & ~noisy(fr) & ~nearSp(fr) & z(n,fr)-base(fr)>=P.spThr & a7(fr)>=P.cs.prom);
    added(n,fr)=true;
end
spAll=double((sp>0 & ~removed) | added);
end

function T=cs_candidates(sub,sp,CS,o,frm_rate)
% envelope transients as in interactive_cs_finder (detect_transient2, spikes within 20 frames
% joined), with the reason each one is or is not a CS
T=struct('on',{},'off',{},'nSp',{},'isiMs',{},'auc',{},'isCS',{},'why',{});
if sum(sp)<2, return; end
[trans,trTrace]=detect_transient2(sub,o.CS_thres,sp,20);
if isempty(trans.amp), return; end
for k=1:numel(trans.spike_number)
    idx=find(trTrace==k);
    T(k).on=idx(1); T(k).off=idx(end);
    T(k).nSp=trans.spike_number(k);
    isiFr=inf;
    if T(k).nSp>1, isiFr=mean(trans.ISI{k}(1:min(o.N_Spike2ISI-1,end))); end
    T(k).isiMs=isiFr/frm_rate*1000;
    T(k).auc=trans.int(k);
    T(k).isCS=any(CS(idx));
    if T(k).isCS, T(k).on=idx(find(CS(idx),1)); T(k).why='';   % span trimmed by the finder
    elseif T(k).nSp<o.N_Spike, T(k).why=sprintf('<%d spikes',o.N_Spike);
    elseif isiFr>=20, T(k).why=sprintf('ISI %.0f ms',T(k).isiMs);
    elseif T(k).auc<=o.Max_AUC, T(k).why=sprintf('AUC %.0f',T(k).auc);
    else, T(k).why='?'; end
end
end

function plot_cs_check(Result,P,figVis)
% complex-spike detection check: traces with CS (orange) and rejected multi-spike
% transients (gray), CS / rejected examples, spikes per CS and first-ISI distributions
f=@(ms) max(1,round(ms*1e-3*Result.frm_rate));
o=Result.CSparams;
[z,sub]=cs_signals(Result,P);
z=z-movmedian(z,500,2);   % display only
t=(1:Result.nFrames)/Result.frm_rate*1000;   % ms
outFile=out_file(Result.fpath,'CScheck.png');
nCS=zeros(size(sub,1),1); Ts=cell(size(sub,1),1);
for n=find(Result.isActive(:))'
    Ts{n}=cs_candidates(sub(n,:),double(Result.spike(n,:)),Result.CStrace(n,:),o,Result.frm_rate);
    if ~isempty(Ts{n}), nCS(n)=sum([Ts{n}.isCS]); end
end
[~,top]=sort(nCS,'descend'); top=top(nCS(top)>0); top=top(1:min(6,end));
fig=figure(302); clf; theme(fig,'light');
set(fig,'Color','w','Position',[50 50 1800 1200],'Visible',figVis);
tiledlayout(5,6,'TileSpacing','compact','Padding','compact');
if isempty(top)
    nexttile([5 6]); text(0.3,0.5,sprintf('%s %s: no complex spikes detected',Result.RecID, ...
        Result.recName),'Interpreter','none'); axis off
    exportgraphics(fig,outFile,'Resolution',110);
    return
end

% traces of the cells with most CS, 3 s around the densest CS period of the top cell
on1=[Ts{top(1)}([Ts{top(1)}.isCS]).on];
c=histcounts(on1,0:f(3000):Result.nFrames); [~,b]=max(c);
w=max(1,(b-1)*f(3000)-f(500)):min(Result.nFrames,b*f(3000)+f(500));
nexttile([2 6]); hold on
for q=1:numel(top)
    n=top(q); y0=(q-1)*25;
    for T=Ts{n}
        if T.off<w(1) || T.on>w(end), continue; end
        col=[1 0.6 0.2]; if ~T.isCS, col=[0.85 0.85 0.85]; end
        if ~T.isCS && T.nSp<2, continue; end
        patch(t([T.on T.off T.off T.on]),y0+[-8 -8 14 14],col,'EdgeColor','none','FaceAlpha',0.5);
    end
    plot(t(w),z(n,w)+y0,'k','LineWidth',0.6);
    plot(t(w),sub(n,w)+y0,'b','LineWidth',1);
    s=w(Result.spike(n,w)); plot(t(s),z(n,s)+y0,'r.','MarkerSize',8);
    s=w(Result.spikeAdded(n,w)); plot(t(s),z(n,s)+y0,'mo','MarkerSize',4);
    s=w(Result.spikeRemoved(n,w)); plot(t(s),z(n,s)+y0,'x','Color',[0.5 0.5 0.5],'MarkerSize',4);
    text(t(w(1)),y0+10,sprintf('cell %d  (%d CS)',n,nCS(n)),'FontSize',8);
end
axis tight; xlabel('time (ms)'); set(gca,'ytick',[]);
title(sprintf(['%s %s | interactive_cs_finder: thres [%g %g], >=%d spikes, ISI < 20 fr, AUC > %g | ' ...
    'spikes on the plateau: %g above pre-plateau baseline, 7-pt amplitude >= %g (magenta: added, gray x: removed ' ...
    'on noisy plateaus) | black: z, blue: envelope, red: spikes | orange: CS, gray: rejected'],Result.RecID,Result.recName, ...
    o.CS_thres,o.N_Spike,o.Max_AUC,P.spThr,P.cs.prom),'Interpreter','none','FontSize',9);

% CS examples and rejected transients (>=2 spikes)
allT=[]; for n=top', for T=Ts{n}, T.cell=n; allT=[allT T]; end, end %#ok<AGROW>
cs=allT([allT.isCS]); rej=allT(~[allT.isCS] & [allT.nSp]>=2);
rng(2); cs=cs(randperm(numel(cs),min(12,numel(cs)))); rej=rej(randperm(numel(rej),min(4,numel(rej))));
ex=[num2cell(cs) num2cell(rej)];
for e=1:16
    nexttile(12+e);
    if e>numel(ex), axis off; continue; end
    T=ex{e}; n=T.cell;
    w=max(1,T.on-f(30)):min(Result.nFrames,T.off+f(60));
    col=[1 0.6 0.2]; if ~T.isCS, col=[0.8 0.8 0.8]; end
    patch(t([T.on T.off T.off T.on])-t(T.on),[-6 -6 20 20],col,'EdgeColor','none','FaceAlpha',0.5); hold on
    plot(t(w)-t(T.on),z(n,w),'k'); plot(t(w)-t(T.on),sub(n,w),'b','LineWidth',1);
    yline(o.CS_thres,':','Color',[0 0 1]);
    s=w(Result.spike(n,w)); plot(t(s)-t(T.on),z(n,s),'r.','MarkerSize',10);
    s=w(Result.spikeAdded(n,w)); plot(t(s)-t(T.on),z(n,s),'mo','MarkerSize',5);
    s=w(Result.spikeRemoved(n,w)); plot(t(s)-t(T.on),z(n,s),'x','Color',[0.5 0.5 0.5],'MarkerSize',5);
    axis tight; ylim([-6 20]);
    if T.isCS
        title(sprintf('CS  cell %d: %d sp, ISI %.0f ms, AUC %.0f',n,T.nSp,T.isiMs,T.auc),'FontSize',7);
    else
        title(sprintf('rejected cell %d: %d sp, %s',n,T.nSp,T.why),'FontSize',7,'Color',[0.5 0 0]);
    end
end

allA=[]; for n=find(Result.isActive(:))', allA=[allA Ts{n}]; end %#ok<AGROW>
nexttile(29);
histogram([allA([allA.isCS]).nSp],2.5:1:12.5,'FaceColor',[1 0.6 0.2]);
xlabel('spikes per CS'); ylabel('# CS'); title(sprintf('%d CS in %d cells',sum([allA.isCS]),sum(nCS>0)),'FontSize',8);
nexttile(30); hold on
m2=[allA.nSp]>=2;
histogram([allA(m2 & [allA.isCS]).isiMs],0:2:60,'FaceColor',[1 0.6 0.2]);
histogram([allA(m2 & ~[allA.isCS]).isiMs],0:2:60,'FaceColor',[0.6 0.6 0.6]);
xline(20/Result.frm_rate*1000,'r--');
xlabel(sprintf('mean first-%d ISI (ms)',o.N_Spike2ISI-1)); legend({'CS','rejected'},'FontSize',7);
exportgraphics(fig,outFile,'Resolution',110);
end

function M=place_tuning(spike,Arduino,valid,frm_rate,nBin)
pos=Arduino(:,2); run=Arduino(:,4)>0 & valid(:);
edges=linspace(min(pos),max(pos)+eps,nBin+1);
pb=discretize(pos,edges);
use=run & ~isnan(pb);
occ=accumarray(pb(use),1,[nBin 1])'/frm_rate;
M=zeros(size(spike,1),nBin);
for n=1:size(spike,1)
    M(n,:)=accumarray(pb(use),double(spike(n,use))',[nBin 1])'./occ;
end
M(:,occ<0.1)=NaN;
end

function plot_recording(Result,P,figVis)
nCell=size(Result.ftprnt,3);
act=find(Result.isActive(:));
col=repmat([0.6 0.6 0.6],nCell,1); col(act,:)=repmat([0 0.45 0.74],numel(act),1);
order=[act; setdiff((1:nCell)',act)];
if nCell>60 && ~isempty(act), order=act; end   % crowded FOVs: show active cells only
t=(1:Result.nFrames)/Result.frm_rate;

fig=figure(300); clf; theme(fig,'light');
set(fig,'Color','w','Position',[50 50 1600 1000],'Visible',figVis);
tiledlayout(12,6,'TileSpacing','compact','Padding','compact');
ax1=nexttile([4 3]);
fp=Result.ftprnt>0.3*max(Result.ftprnt,[],[1 2]);
show_footprnt_contour(fp,Result.ref_im,col);
hold on
cmapP=[0 0.45 0.74; 0.85 0.33 0.1; 0.93 0.69 0.13; 0.49 0.18 0.56];
for j=1:numel(Result.DMDPatt)
    if isempty(Result.DMDPatt{j}), continue; end
    plot(Result.DMDPatt{j}(:,1),Result.DMDPatt{j}(:,2),'--','Color',cmapP(mod(j-1,4)+1,:),'LineWidth',1);
end
title(sprintf('%d footprints, %d active (blue)',nCell,numel(act)));

nexttile([4 3]);
if ~isempty(Result.placeMap)
    pm=Result.placeMap(Result.isActive,:);
    pm=pm./max(pm,[],2);
    [~,pk]=max(pm,[],2); [~,o]=sort(pk);
    imagesc(pm(o,:)); colormap(gca,turbo); xlabel('position bin'); ylabel('active cells');
    title('Place tuning (running, sorted by peak)');
elseif ~isempty(act)
    histogram(Result.spRate(act),'FaceColor',[0 0.45 0.74]);
    xlabel('firing rate (Hz)'); ylabel('# active cells');
    title(sprintf('Firing rate, median %.1f Hz',median(Result.spRate(act))));
else
    axis off
end

ax2=nexttile([6 6]);
tr=min(max(double(Result.zTr(order,:)),-6),20);
ds=max(1,floor(Result.nFrames/30000));   % max-decimate long recordings for display
if ds>1
    nb=floor(Result.nFrames/ds);
    trD=squeeze(max(reshape(tr(:,1:nb*ds),numel(order),ds,nb),[],2));
    tD=t(ds:ds:nb*ds);
else
    trD=tr; tD=t;
end
cs=Result.CStrace(order,:);
if ds>1, csD=squeeze(any(reshape(cs(:,1:nb*ds),numel(order),ds,nb),2)); else, csD=cs; end
scale=12;
hold on
for n=1:numel(order)
    plot(tD,trD(n,:)+n*scale,'Color',col(order(n),:),'LineWidth',0.5);
    y=trD(n,:); y(~csD(n,:))=NaN;
    plot(tD,y+n*scale,'r','LineWidth',0.9);   % complex spikes
    s=find(Result.spike(order(n),:) & ~cs(n,:));
    plot(t(s),tr(n,s)+n*scale,'k.','MarkerSize',4);
end
set(gca,'ytick',(1:numel(order))*scale,'yticklabel',order,'FontSize',5);
axis tight; ylabel('cell');
title(sprintf('%s  %s  |  %s: %s  |  %s  |  black dots: spikes, red: complex spikes (%d)', ...
    Result.RecID,Result.recName,Result.protocol,Result.protocolDetail,Result.state,sum(Result.nCS)), ...
    'Interpreter','none','FontSize',9);

ax3=nexttile([2 6]); hold on
if ~isempty(Result.Arduino)
    yyaxis left; plot(t,Result.Arduino(:,2),'k'); ylabel('position');
    yyaxis right; plot(t,Result.Blue,'b'); ylabel('blue (V)');
elseif ~isempty(Result.BluePatt)
    for j=1:size(Result.BluePatt,1)
        plot(t,Result.BluePatt(j,:),'Color',cmapP(mod(j-1,4)+1,:));
    end
    ylabel('blue per DMD pattern (V)');
else
    plot(t,Result.Blue,'b'); ylabel('blue (V)');
end
if any(Result.excluded)   % laser off, motion transients, global artifacts
    yl=ylim; area(t,Result.excluded*yl(2),'FaceColor',[1 0.6 0.6],'EdgeColor','none','FaceAlpha',0.4);
end
xlabel('time (s)'); axis tight
linkaxes([ax2 ax3],'x');
exportgraphics(fig,out_file(Result.fpath,'traces.png'),'Resolution',150);
end
