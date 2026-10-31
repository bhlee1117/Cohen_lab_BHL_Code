function [cfg, R] = flow_gui(data, cfg, R)
% [cfg, R] = flow_gui(data)
% [cfg, R] = flow_gui(data, cfg)
% [cfg, R] = flow_gui(data, cfg, R)      % resume from an earlier session
%
% Interactive front end for the vascular-flow pipeline: the six stages across the top,
% the selected stage's parameters in an editable table, and its result in the preview.
% Adjust a value, run that stage, look, adjust again.
%
% Blocks until you press Done (or close the window) and then returns the parameters and
% the results, so it can be used inline in a script:
%
%     [cfg, R] = flow_gui(dataFlow);         % tune interactively, then carry on
%     save(fullfile(outDir,'flow_cfg.mat'), 'cfg');
%
% INPUT
%   data  the recording, as built for flow_pipeline:
%     .mov_res [nrow x ncol x nT] required. .fs .pxUm .badfr .spike .ftprnt .bvBW
%     .neuronMap .voltTr .tag       (missing optional fields are filled with defaults)
%   cfg   parameter set; [] or omitted = flow_defaults()
%   R     previous results to resume from ([] = start fresh)
%
% OUTPUT
%   cfg   the parameters as left in the table -- feed straight back in next time
%   R     R.<stage> for every stage that ran, plus R.stamp.<stage>
%
% HOW STALENESS IS SHOWN. Each step button is coloured by state:
%   grey    not run           green  up to date
%   amber   STALE: a parameter it depends on changed since it ran
%   red     the stage function does not exist yet (not ported)
% Amber is the point of the whole thing. In the script version nothing tracked which
% cells a parameter change invalidated, so the choice was to re-run everything or to
% re-run one cell and hope. Here the stamp chain answers it, and running a stage clears
% the amber on that stage and sets it on everything downstream.
%
% THE PARAMETER TABLE is generated from the cfg struct itself, not hand-written, so a
% knob added to flow_defaults appears here with no change to this file. Values are typed
% as MATLAB literals: 40, [200 220], 'dog', true, {'full','fast'}.
%
% 2026.08.28, Byung Hun Lee

if nargin < 2 || isempty(cfg), cfg = flow_defaults(); end
if nargin < 3, R = struct; end
assert(isstruct(data) && isfield(data,'mov_res'), 'data.mov_res is required.');

% Validate and fill in once, up front: a mask that does not match the movie should be
% reported by name here, not as a size error inside a stage. The same goes for the config
% -- a cfg saved before a knob existed used to fail as "Unrecognized field name" deep
% inside whichever stage first read it, which says nothing about the actual problem.
data = flow_checkData(data);
cfg  = flow_checkCfg(cfg);
sz = [size(data.mov_res,1) size(data.mov_res,2)];  nT = size(data.mov_res,3);

S.stages = {'bandpass','track','segments','regressors','apply','compare'};
S.label  = {'1 bandpass','2 track','3 segments','4 regressors','5 apply','6 compare'};
S.data = data;  S.cfg = cfg;  S.R = R;  S.sel = 1;

fh = figure('Name','FLOW pipeline','NumberTitle','off','Color',[0.94 0.94 0.94], ...
            'Units','normalized','Position',[0.08 0.10 0.80 0.78], ...
            'CloseRequestFcn',@(~,~) finish());
flow_theme(fh);          % light, and every figure a stage opens gets it too

% ---- step bar -----------------------------------------------------------------------
nS = numel(S.stages);
for k = 1:nS
    S.btn(k) = uicontrol(fh,'Style','pushbutton','String',S.label{k}, ...
        'Units','normalized','Position',[0.015+(k-1)*0.1625 0.930 0.152 0.055], ...
        'FontSize',10,'Callback',@(~,~) selectStage(k));
end

% ---- actions ------------------------------------------------------------------------
uicontrol(fh,'Style','pushbutton','String','Run this step','Units','normalized', ...
    'Position',[0.015 0.872 0.115 0.048],'FontSize',10,'Callback',@(~,~) runTo(false));
uicontrol(fh,'Style','pushbutton','String','Run up to here','Units','normalized', ...
    'Position',[0.135 0.872 0.115 0.048],'FontSize',10,'Callback',@(~,~) runTo(true));
uicontrol(fh,'Style','pushbutton','String','Reset params','Units','normalized', ...
    'Position',[0.255 0.872 0.100 0.048],'Callback',@(~,~) resetParams());
uicontrol(fh,'Style','pushbutton','String','Done','Units','normalized', ...
    'Position',[0.885 0.872 0.100 0.048],'FontSize',10,'FontWeight','bold', ...
    'Callback',@(~,~) finish());
S.status = uicontrol(fh,'Style','text','String','','Units','normalized', ...
    'Position',[0.365 0.875 0.50 0.042],'HorizontalAlignment','left', ...
    'BackgroundColor',get(fh,'Color'),'FontSize',9);

% ---- LEFT: what this step does, its parameters, then the log -------------------------
S.info = flow_paramInfo();
S.desc = uicontrol(fh,'Style','text','String','','Units','normalized', ...
    'Position',[0.015 0.745 0.42 0.112],'HorizontalAlignment','left', ...
    'FontSize',9.5,'BackgroundColor',[1 1 1]);
% The note column is read out of flow_defaults.m, so a knob's explanation here is the
% same text as the comment beside its value -- they cannot drift apart.
S.tbl = uitable(fh,'Units','normalized','Position',[0.015 0.360 0.42 0.375], ...
    'ColumnName',{'parameter','value','what it does'},'ColumnEditable',[false true false], ...
    'ColumnWidth',{125 90 'auto'},'RowName',[], ...
    'CellEditCallback',@(src,e) onEdit(e));
S.logAx = axes(fh,'Units','normalized','Position',[0.015 0.025 0.42 0.315]);
S.logLines = flow_log(S.logAx, [], 'init');

% ---- RIGHT: the stage's own views ----------------------------------------------------
S.panel = uipanel(fh,'Units','normalized','Position',[0.450 0.025 0.535 0.832], ...
    'BackgroundColor',[1 1 1],'BorderType','line');

guidata(fh, S);
refresh_();
say(sprintf(['movie %dx%dx%d, fs %g Hz, %g um/px\n' ...
             'pick a step, edit its parameters, press Run this step.'], ...
             sz(1), sz(2), nT, data.fs, data.pxUm));
uiwait(fh);

if isvalid(fh)                      % Done pressed: read back, then close for real
    S = guidata(fh);  cfg = S.cfg;  R = S.R;
    delete(fh);
else
    cfg = S.cfg;  R = S.R;
end

% =====================================================================================
    function selectStage(k)
        S = guidata(fh);  S.sel = k;  guidata(fh, S);  refresh_();
    end

    function onEdit(e)
        S = guidata(fh);
        st = S.stages{S.sel};
        nm = S.tbl.Data{e.Indices(1), 1};
        [v, ok, msg] = flow_str2val(e.NewData);
        if ~ok
            say(sprintf('could not parse "%s" for %s: %s', e.NewData, nm, msg), true);
            S.tbl.Data{e.Indices(1),2} = flow_val2str(S.cfg.(st).(nm));   % put the old one back
            return
        end
        S.cfg.(st).(nm) = v;
        guidata(fh, S);
        say(sprintf('%s.%s = %s', st, nm, flow_val2str(v)));
        refresh_();                 % staleness may have changed downstream
    end

    function resetParams()
        S = guidata(fh);
        dflt = flow_defaults();
        st = S.stages{S.sel};
        S.cfg.(st) = dflt.(st);
        guidata(fh, S);
        say(sprintf('%s parameters reset to defaults', st));
        refresh_();
    end

    function runTo(chain)
        S = guidata(fh);
        st = S.stages{S.sel};
        miss = missingStages(S, S.sel, chain);
        if ~isempty(miss)
            say(sprintf(['cannot run: %s not ported yet.\n' ...
                'The stage functions still to be written are listed in flow_pipeline.'], ...
                strjoin(miss, ', ')), true);
            return
        end
        set(S.status,'String',sprintf('running %s ...', st));  drawnow;
        t0 = tic;
        try
            out = evalc('S.R = flow_pipeline(S.data, S.cfg, S.R, st);');
            guidata(fh, S);
            say(sprintf('%s\n-- %s finished in %.1f s --', strtrim(out), st, toc(t0)));
        catch ME
            loc = '';
            if ~isempty(ME.stack), loc = sprintf(' (%s line %d)', ME.stack(1).name, ME.stack(1).line); end
            say(sprintf('%s FAILED: %s%s', st, ME.message, loc), true);
        end
        set(S.status,'String','');
        refresh_();
    end

    function finish()
        if isvalid(fh), uiresume(fh); end
    end

    function say(txt, isErr)
        if nargin < 2, isErr = false; end
        S = guidata(fh);
        S.logLines = flow_log(S.logAx, S.logLines, 'add', txt, isErr);
        guidata(fh, S);
    end

    function refresh_()
        S = guidata(fh);
        st = S.stages{S.sel};

        % --- parameter table for the selected stage ---
        set(S.desc, 'String', sprintf('STEP %d  %s\n\n%s', S.sel, upper(st), ...
            getfld(S.info.stage, st, '')));
        p = S.cfg.(st);  f2 = fieldnames(p);
        D = cell(numel(f2), 3);
        for q = 1:numel(f2)
            note = '';
            if isfield(S.info.param, st), note = getfld(S.info.param.(st), f2{q}, ''); end
            % say so when a knob is inert given the other settings, instead of leaving
            % the user to discover it by changing it and seeing nothing happen
            if strcmp(st,'bandpass') && strcmp(f2{q},'order') && ~strcmpi(p.method,'fourier')
                note = ['[INACTIVE while method = ''' p.method '''] ' note];
            end
            if strcmp(st,'bandpass') && startsWith(f2{q},'ring') && ~strcmpi(p.voltRef,'ring')
                note = ['[INACTIVE while voltRef = ''' p.voltRef '''] ' note];
            end
            D{q,1} = f2{q};  D{q,2} = flow_val2str(p.(f2{q}));  D{q,3} = note;
        end
        set(S.tbl, 'Data', D);

        % --- colour every step button by its state ---
        want = wantStamps(S.cfg, S.stages);
        for q = 1:numel(S.stages)
            sq = S.stages{q};
            if isempty(which(['flow_' sq]))
                col = [0.92 0.62 0.62];  tag = 'not ported';
            elseif ~isfield(S.R, sq)
                col = [0.86 0.86 0.86];  tag = 'not run';
            elseif isfield(S.R,'stamp') && isfield(S.R.stamp, sq) && ...
                   strcmp(S.R.stamp.(sq), want.(sq))
                col = [0.72 0.90 0.72];  tag = 'up to date';
            else
                col = [0.98 0.85 0.60];  tag = 'STALE';
            end
            set(S.btn(q), 'BackgroundColor', col, ...
                'String', sprintf('%s\n%s', S.label{q}, tag));
        end
        set(S.btn(S.sel), 'FontWeight','bold');
        set(S.btn(setdiff(1:numel(S.stages), S.sel)), 'FontWeight','normal');

        % --- preview: the stage lays out its own subplots inside the panel ---
        flow_preview(S.panel, st, S.R, S.cfg, S.data);
        guidata(fh, S);
    end
end

% =========================================================================
function v = getfld(s, f, dflt)
if isstruct(s) && isfield(s, f), v = s.(f); else, v = dflt; end
end

function miss = missingStages(S, upto, chain)
% which of the stages needed for this run do not exist as functions yet
if chain, idx = 1:upto; else, idx = 1:upto; end   % a single step still needs its upstream
miss = {};
for k = idx
    if isempty(which(['flow_' S.stages{k}])), miss{end+1} = S.stages{k}; end %#ok<AGROW>
end
end

function want = wantStamps(cfg, stages)
% the stamp each stage SHOULD have, chained exactly as flow_pipeline does
prev = flow_stamp(cfg.bandpass, cfg.common);
want.(stages{1}) = prev;
for k = 2:numel(stages)
    prev = flow_stamp(cfg.(stages{k}), prev);
    want.(stages{k}) = prev;
end
end
