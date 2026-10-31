function dataOrganizer(mode)
% DATAORGANIZER  GUI to collect data files scattered across several disks.
%
%   dataOrganizer              opens the GUI
%   dataOrganizer('selftest')  runs a self-test on a temporary folder tree
%
% Workflow:
%   1. Add one or more source folders (they can live on different disks).
%   2. Set the file names to look for (default: frame1.bin, frame2.bin,
%      output_data.mat). Wildcards are allowed, e.g. frame*.bin or *.mat.
%   3. Pick a destination folder.
%   4. Press SCAN to list every match found in the sources and all of their
%      subfolders, then review / untick rows in the table.
%   5. Press RUN to move (or copy) the ticked files.
%
% Safety notes:
%   * Nothing is overwritten unless you pick "Overwrite" as the conflict policy.
%   * A move across disks is copy -> verify size -> delete source.  A move
%     inside one disk is a fast rename.
%   * Every real run writes a manifest (.json) into the destination folder,
%     and "Undo from manifest..." puts the files back where they came from.
%   * "Dry run" logs everything it would do without touching any file.

if nargin > 0 && strcmpi(char(mode),'selftest')
    selfTest;
    return;
end

%% ------------------------------------------------------------------ state
S = struct();
S.roots        = {};          % cellstr of source folders
S.items        = emptyItems;  % struct array of found files
S.lastManifest = '';

%% --------------------------------------------------------------------- UI
fig = uifigure('Name','Data Organizer','Position',[100 100 1200 780]);
fig.CloseRequestFcn = @onClose;

outer = uigridlayout(fig,[3 1]);
outer.RowHeight   = {215,'1x',170};
outer.ColumnWidth = {'1x'};

top = uigridlayout(outer,[1 3]);
top.Layout.Row    = 1;
top.Layout.Column = 1;
top.ColumnWidth   = {'1.1x','1x','1.15x'};
top.Padding       = [0 0 0 0];

% ---------------------------------------------------- 1. sources
pSrc = uipanel(top,'Title','1. Source folders  (add one per disk / dataset)');
pSrc.Layout.Row = 1; pSrc.Layout.Column = 1;
gSrc = uigridlayout(pSrc,[2 4]);
gSrc.RowHeight   = {'1x',26};
gSrc.ColumnWidth = {'1x','fit','fit','fit'};

lstRoots = uilistbox(gSrc,'Items',{},'Multiselect','on');
lstRoots.Layout.Row = 1; lstRoots.Layout.Column = [1 4];

b = uibutton(gSrc,'Text','Add folder...','ButtonPushedFcn',@(~,~)addRoot);
b.Layout.Row = 2; b.Layout.Column = 1;
b = uibutton(gSrc,'Text','Remove','ButtonPushedFcn',@(~,~)removeRoot);
b.Layout.Row = 2; b.Layout.Column = 2;
b = uibutton(gSrc,'Text','Clear','ButtonPushedFcn',@(~,~)clearRoots);
b.Layout.Row = 2; b.Layout.Column = 3;
b = uibutton(gSrc,'Text','Show disks','ButtonPushedFcn',@(~,~)showDisks);
b.Layout.Row = 2; b.Layout.Column = 4;

% ---------------------------------------------------- 2. what to move
pWhat = uipanel(top,'Title','2. Files to look for');
pWhat.Layout.Row = 1; pWhat.Layout.Column = 2;
gWhat = uigridlayout(pWhat,[4 1]);
gWhat.RowHeight   = {'fit','1x','fit','fit'};
gWhat.ColumnWidth = {'1x'};

lb = uilabel(gWhat,'Text','Names / wildcards (comma or new line separated):');
lb.Layout.Row = 1; lb.Layout.Column = 1;

txtPat = uitextarea(gWhat,'Value',{'frame1.bin, frame2.bin, output_data.mat'});
txtPat.Layout.Row = 2; txtPat.Layout.Column = 1;

chkHidden = uicheckbox(gWhat,'Text','Skip hidden / system folders','Value',true);
chkHidden.Layout.Row = 3; chkHidden.Layout.Column = 1;

chkNoLinks = uicheckbox(gWhat,'Text','Ignore symlinks / aliases','Value',true);
chkNoLinks.Layout.Row = 4; chkNoLinks.Layout.Column = 1;

% ---------------------------------------------------- 3. destination
pDst = uipanel(top,'Title','3. Destination and options');
pDst.Layout.Row = 1; pDst.Layout.Column = 3;
gDst = uigridlayout(pDst,[7 3]);
gDst.RowHeight   = repmat({22},1,7);
gDst.ColumnWidth = {80,'1x',95};
gDst.RowSpacing  = 4;

lb = uilabel(gDst,'Text','Folder:');                lb.Layout.Row = 1; lb.Layout.Column = 1;
txtDst = uieditfield(gDst,'text','Value','');       txtDst.Layout.Row = 1; txtDst.Layout.Column = 2;
b = uibutton(gDst,'Text','Browse...','ButtonPushedFcn',@(~,~)browseDest);
b.Layout.Row = 1; b.Layout.Column = 3;

lb = uilabel(gDst,'Text','Operation:');             lb.Layout.Row = 2; lb.Layout.Column = 1;
ddOp = uidropdown(gDst,'Items',{'Move','Copy'},'Value','Move');
ddOp.Layout.Row = 2; ddOp.Layout.Column = [2 3];

lb = uilabel(gDst,'Text','Layout:');                lb.Layout.Row = 3; lb.Layout.Column = 1;
ddLayout = uidropdown(gDst, ...
    'Items',{'Mirror source subfolders','Flat, prefix parent folder names','Flat, keep original names'}, ...
    'Value','Mirror source subfolders','ValueChangedFcn',@(~,~)syncLayoutOpts);
ddLayout.Layout.Row = 3; ddLayout.Layout.Column = [2 3];

lb = uilabel(gDst,'Text','Prefix depth:');          lb.Layout.Row = 4; lb.Layout.Column = 1;
spnDepth = uispinner(gDst,'Limits',[1 6],'Value',1,'Step',1,'Enable','off');
spnDepth.Layout.Row = 4; spnDepth.Layout.Column = 2;
lblDepth = uilabel(gDst,'Text','parent folder(s)','Enable','off');
lblDepth.Layout.Row = 4; lblDepth.Layout.Column = 3;

lb = uilabel(gDst,'Text','On conflict:');           lb.Layout.Row = 5; lb.Layout.Column = 1;
ddConf = uidropdown(gDst,'Items',{'Auto-rename (_1, _2, ...)','Skip file','Overwrite'}, ...
    'Value','Auto-rename (_1, _2, ...)','ValueChangedFcn',@(~,~)replan);
ddConf.Layout.Row = 5; ddConf.Layout.Column = [2 3];

chkRootName = uicheckbox(gDst,'Text','Include source folder name as top level', ...
    'Value',true,'ValueChangedFcn',@(~,~)replan);
chkRootName.Layout.Row = 6; chkRootName.Layout.Column = [1 3];

chkDry = uicheckbox(gDst,'Text','Dry run (log only, change nothing)','Value',false);
chkDry.Layout.Row = 7; chkDry.Layout.Column = [1 3];

% ---------------------------------------------------- 4. preview
pTab = uipanel(outer,'Title','4. Preview');
pTab.Layout.Row = 2; pTab.Layout.Column = 1;
gTab = uigridlayout(pTab,[2 1]);
gTab.RowHeight   = {28,'1x'};
gTab.ColumnWidth = {'1x'};

bar = uigridlayout(gTab,[1 9]);
bar.Layout.Row = 1; bar.Layout.Column = 1;
bar.ColumnWidth = {110,70,70,80,'1x',150,110,110,180};
bar.Padding     = [0 0 0 0];

uibutton(bar,'Text','SCAN','FontWeight','bold','ButtonPushedFcn',@(~,~)doScan);
uibutton(bar,'Text','All','ButtonPushedFcn',@(~,~)setAll(true));
uibutton(bar,'Text','None','ButtonPushedFcn',@(~,~)setAll(false));
uibutton(bar,'Text','Invert','ButtonPushedFcn',@(~,~)invertAll);
lblCount = uilabel(bar,'Text','No scan yet.');
uibutton(bar,'Text','Reveal in Finder','ButtonPushedFcn',@(~,~)revealSelected);
uibutton(bar,'Text','Export plan','ButtonPushedFcn',@(~,~)exportCsv);
uibutton(bar,'Text','RUN','FontWeight','bold','BackgroundColor',[0.30 0.62 0.36], ...
    'FontColor','white','ButtonPushedFcn',@(~,~)doRun);
uibutton(bar,'Text','Undo from manifest...','ButtonPushedFcn',@(~,~)doUndo);

tbl = uitable(gTab);
tbl.Layout.Row = 2; tbl.Layout.Column = 1;
tbl.ColumnName       = {'Move','File','Size','Source folder','Destination','Status'};
tbl.ColumnEditable   = [true false false false false false];
tbl.ColumnWidth      = {45,190,80,'auto','auto',160};
tbl.CellEditCallback = @onCellEdit;
tbl.Data             = cell(0,6);

% ---------------------------------------------------- log
pLog = uipanel(outer,'Title','Log');
pLog.Layout.Row = 3; pLog.Layout.Column = 1;
gLog = uigridlayout(pLog,[1 1]);
txtLog = uitextarea(gLog,'Value',{''},'Editable','off','FontName','Menlo');
txtLog.Layout.Row = 1; txtLog.Layout.Column = 1;

loadPrefs;
syncLayoutOpts;
log_('Ready. Add source folders, set the destination, then press SCAN.');

%% ------------------------------------------------------------- source list
    function addRoot
        d = uigetdir_(lastDirGuess,'Select a source folder');
        figure(fig);
        if isempty(d), return; end
        if any(strcmp(S.roots,d))
            log_('Already in the list: %s', d); return;
        end
        S.roots{end+1} = d;
        refreshRoots;
        log_('Added source: %s', d);
    end

    function removeRoot
        sel = lstRoots.Value;
        if isempty(sel), return; end
        if ~iscell(sel), sel = {sel}; end
        S.roots(ismember(S.roots,sel)) = [];
        refreshRoots;
    end

    function clearRoots
        S.roots = {};
        refreshRoots;
    end

    function refreshRoots
        lstRoots.Items = S.roots;
    end

    function showDisks
        if ismac
            d = dir('/Volumes');
            names = {d([d.isdir]).name};
            names = names(~ismember(names,{'.','..'}));
            log_('Mounted volumes: %s', strjoin(names,'  |  '));
        elseif isunix
            log_('Look under /media or /mnt for mounted disks.');
        else
            log_('Use "Add folder..." to browse the available drives.');
        end
    end

    function d = lastDirGuess
        if ~isempty(S.roots),          d = S.roots{end};
        elseif ~isempty(txtDst.Value), d = txtDst.Value;
        else,                          d = pwd;
        end
    end

    function browseDest
        d = uigetdir_(txtDst.Value,'Select the destination folder');
        figure(fig);
        if ~isempty(d)
            txtDst.Value = d;
            replan;
        end
    end

    function syncLayoutOpts
        isPrefix = strcmp(ddLayout.Value,'Flat, prefix parent folder names');
        isMirror = strcmp(ddLayout.Value,'Mirror source subfolders');
        onoff = {'off','on'};
        spnDepth.Enable    = onoff{isPrefix+1};
        lblDepth.Enable    = onoff{isPrefix+1};
        chkRootName.Enable = onoff{isMirror+1};
        replan;
    end

    function replan
        if isempty(S.items), return; end
        planDestinations;
        refreshTable;
    end

%% -------------------------------------------------------------------- scan
    function doScan
        dst  = strtrim(txtDst.Value);
        pats = parsePatterns(txtPat.Value);
        if isempty(S.roots)
            uialert(fig,'Add at least one source folder.','No sources'); return;
        end
        if isempty(pats)
            uialert(fig,'Enter at least one file name or wildcard.','No patterns'); return;
        end

        roots = dropNestedRoots(S.roots);
        S.items = emptyItems;
        refreshTable;

        dlg = uiprogressdlg(fig,'Title','Scanning','Message','Searching subfolders...', ...
                            'Indeterminate',true,'Cancelable',true);
        guard = onCleanup(@()closeIfValid(dlg));

        found = emptyItems;
        seen  = containers.Map('KeyType','char','ValueType','logical');
        t0 = tic;

        for ir = 1:numel(roots)
            root = roots{ir};
            if ~isfolder(root)
                log_('SKIP (not reachable): %s', root);
                continue;
            end
            for ip = 1:numel(pats)
                if dlg.CancelRequested, break; end
                dlg.Message = sprintf('%s\nlooking for %s', root, pats{ip});
                drawnow limitrate;
                try
                    d = dir(fullfile(root,'**',pats{ip}));
                catch ME
                    log_('Scan error in %s: %s', root, ME.message);
                    continue;
                end
                if isempty(d), continue; end
                d = d(~[d.isdir]);
                for k = 1:numel(d)
                    full = fullfile(d(k).folder, d(k).name);
                    if isKey(seen,full), continue; end
                    if chkHidden.Value && isHiddenPath(d(k).folder, root, d(k).name), continue; end
                    if chkNoLinks.Value && isLink(full), continue; end
                    if ~isempty(dst) && isUnder(full,dst), continue; end  % already at destination
                    seen(full) = true;
                    found(end+1) = struct('src',full,'name',d(k).name,'bytes',d(k).bytes, ...
                                          'root',root,'dst','','status','','include',true); %#ok<AGROW>
                end
            end
            if dlg.CancelRequested
                log_('Scan cancelled after %d source(s).', ir-1);
                break;
            end
        end

        S.items = found;
        planDestinations;
        refreshTable;
        log_('Scan finished in %s: %d file(s), %s total.', ...
             humanTime(toc(t0)), numel(S.items), humanBytes(totalBytes(S.items)));
        if isempty(S.items)
            log_('Nothing matched. Check the names/wildcards and that the disks are mounted.');
        end
        savePrefs;
    end

%% ---------------------------------------------------- destination planning
    function planDestinations
        if isempty(S.items), return; end
        dst = strtrim(txtDst.Value);
        if isempty(dst)
            for i = 1:numel(S.items)
                S.items(i).dst    = '';
                S.items(i).status = 'No destination set';
            end
            return;
        end

        mode   = ddLayout.Value;
        policy = ddConf.Value;
        taken  = containers.Map('KeyType','char','ValueType','logical');

        for i = 1:numel(S.items)
            it = S.items(i);
            switch mode
                case 'Mirror source subfolders'
                    rel = relPath(it.src, it.root);      % e.g. expA/run3/frame1.bin
                    if chkRootName.Value
                        target = fullfile(dst, safeName(baseName(it.root)), rel);
                    else
                        target = fullfile(dst, rel);
                    end
                case 'Flat, prefix parent folder names'
                    target = fullfile(dst, [prefixOf(it.src, spnDepth.Value) it.name]);
                otherwise
                    target = fullfile(dst, it.name);
            end

            status = 'Ready';
            if isKey(taken, lower(target)) || isfile(target)
                switch policy
                    case 'Auto-rename (_1, _2, ...)'
                        target = uniqueTarget(target, taken);
                        status = 'Ready (renamed)';
                    case 'Skip file'
                        status = 'Skip (exists)';
                    otherwise
                        status = 'Ready (overwrite)';
                end
            end
            taken(lower(target)) = true;
            S.items(i).dst    = target;
            S.items(i).status = status;
        end
    end

%% ------------------------------------------------------------------- table
    function refreshTable
        n = numel(S.items);
        D = cell(n,6);
        for i = 1:n
            it = S.items(i);
            D{i,1} = logical(it.include);
            D{i,2} = it.name;
            D{i,3} = humanBytes(it.bytes);
            D{i,4} = fileparts(it.src);
            D{i,5} = it.dst;
            D{i,6} = it.status;
        end
        tbl.Data = D;
        if n == 0
            lblCount.Text = 'No files listed.';
        else
            inc = logical([S.items.include]);
            lblCount.Text = sprintf('%d of %d ticked  (%s)', sum(inc), n, ...
                                    humanBytes(totalBytes(S.items(inc))));
        end
    end

    function onCellEdit(~, ev)
        r = ev.Indices(1);
        c = ev.Indices(2);
        if c == 1 && r >= 1 && r <= numel(S.items)
            S.items(r).include = logical(ev.NewData);
            refreshTable;
        end
    end

    function setAll(tf)
        for i = 1:numel(S.items), S.items(i).include = tf; end
        refreshTable;
    end

    function invertAll
        for i = 1:numel(S.items), S.items(i).include = ~S.items(i).include; end
        refreshTable;
    end

    function revealSelected
        sel = tbl.Selection;
        if isempty(sel), log_('Select a row in the table first.'); return; end
        r = sel(1,1);
        if r < 1 || r > numel(S.items), return; end
        p = fileparts(S.items(r).src);
        if ismac
            system(sprintf('open "%s"', p));
        elseif ispc
            winopen(p);
        else
            system(sprintf('xdg-open "%s" &', p));
        end
    end

    function exportCsv
        if isempty(S.items), log_('Nothing to export.'); return; end
        [f,p] = uiputfile('data_organizer_plan.csv','Save the plan as');
        figure(fig);
        if isequal(f,0), return; end
        T = table({S.items.src}',{S.items.dst}',[S.items.bytes]', ...
                  {S.items.status}',[S.items.include]', ...
                  'VariableNames',{'Source','Destination','Bytes','Status','Ticked'});
        writetable(T, fullfile(p,f));
        log_('Plan written to %s', fullfile(p,f));
    end

%% --------------------------------------------------------------------- run
    function doRun
        dst = strtrim(txtDst.Value);
        if isempty(S.items), uialert(fig,'Run SCAN first.','Nothing to do'); return; end
        if isempty(dst),     uialert(fig,'Pick a destination folder.','No destination'); return; end

        planDestinations;          % re-plan in case options changed after the scan
        refreshTable;

        idx = find([S.items.include] & ~strcmp({S.items.status},'Skip (exists)'));
        if isempty(idx), uialert(fig,'No ticked files to process.','Nothing to do'); return; end

        isMove  = strcmp(ddOp.Value,'Move');
        dry     = chkDry.Value;
        totalB  = totalBytes(S.items(idx));
        verb    = ddOp.Value;
        if dry, verb = [verb ' (dry run)']; end

        q = uiconfirm(fig, sprintf('%s %d file(s), %s, into\n%s', ...
                verb, numel(idx), humanBytes(totalB), dst), 'Confirm', ...
                'Options',{'Go ahead','Cancel'},'DefaultOption',2,'CancelOption',2,'Icon','question');
        if ~strcmp(q,'Go ahead'), return; end

        if ~dry && ~isfolder(dst)
            [ok,msg] = mkdir(dst);
            if ~ok
                uialert(fig,sprintf('Cannot create the destination:\n%s',msg),'Error'); return;
            end
        end

        dlg = uiprogressdlg(fig,'Title',verb,'Message','Starting...','Cancelable',true,'Value',0);
        guard = onCleanup(@()closeIfValid(dlg));

        nOk = 0; nFail = 0; doneB = 0;
        manifest = struct('op',{},'src',{},'dst',{},'bytes',{});
        t0 = tic;
        log_('--- %s started: %d file(s), %s ---', upper(verb), numel(idx), humanBytes(totalB));

        for j = 1:numel(idx)
            i  = idx(j);
            it = S.items(i);
            if dlg.CancelRequested
                log_('Cancelled by user after %d of %d file(s).', j-1, numel(idx));
                break;
            end
            dlg.Value   = min(1, doneB/max(totalB,1));
            dlg.Message = sprintf('%d / %d   %s', j, numel(idx), it.name);
            drawnow limitrate;

            if dry
                log_('[dry] %s  ->  %s', it.src, it.dst);
                S.items(i).status = 'Dry run OK';
                nOk = nOk + 1;
                doneB = doneB + it.bytes;
                continue;
            end

            [ok,msg] = transferOne(it, isMove, strcmp(ddConf.Value,'Overwrite'));
            if ok
                S.items(i).status = 'Done';
                nOk = nOk + 1;
                manifest(end+1) = struct('op',lower(ddOp.Value),'src',it.src, ...
                                         'dst',it.dst,'bytes',it.bytes); %#ok<AGROW>
            else
                S.items(i).status = ['FAILED: ' msg];
                nFail = nFail + 1;
                log_('FAILED  %s  :  %s', it.src, msg);
            end
            doneB = doneB + it.bytes;
        end

        closeIfValid(dlg);
        refreshTable;
        nSkip = sum(strcmp({S.items.status},'Skip (exists)'));
        log_('--- finished in %s: %d ok, %d failed, %d skipped (already present) ---', ...
             humanTime(toc(t0)), nOk, nFail, nSkip);

        if ~dry && ~isempty(manifest)
            mf = fullfile(dst, sprintf('dataOrganizer_manifest_%s.json', stamp('yyyyMMdd_HHmmss')));
            try
                fid = fopen(mf,'w');
                if fid < 0, error('cannot open %s', mf); end
                fprintf(fid,'%s',jsonencode(manifest,'PrettyPrint',true));
                fclose(fid);
                S.lastManifest = mf;
                log_('Manifest written: %s', mf);
            catch ME
                log_('Could not write the manifest: %s', ME.message);
            end
        end
    end

%% -------------------------------------------------------------------- undo
    function doUndo
        start = S.lastManifest;
        if isempty(start), start = strtrim(txtDst.Value); end
        if isempty(start), start = pwd; end
        [f,p] = uigetfile({'*.json','Manifest (*.json)'},'Pick a manifest to undo', start);
        figure(fig);
        if isequal(f,0), return; end
        mf = fullfile(p,f);
        try
            m = jsondecode(fileread(mf));
        catch ME
            uialert(fig,sprintf('Cannot read the manifest:\n%s',ME.message),'Error'); return;
        end
        if iscell(m), m = [m{:}]; end
        if ~isstruct(m) || ~all(isfield(m,{'op','src','dst'}))
            uialert(fig,'That file is not a Data Organizer manifest.','Error'); return;
        end

        nMove = sum(strcmp({m.op},'move'));
        nCopy = numel(m) - nMove;
        q = uiconfirm(fig, sprintf(['Undo %d entry(ies) from\n%s\n\n' ...
                '%d moved file(s) would go back to their original folder.\n' ...
                '%d copied file(s) would be deleted at the destination.'], ...
                numel(m), mf, nMove, nCopy), 'Undo', ...
                'Options',{'Undo moves only','Undo moves + delete copies','Cancel'}, ...
                'DefaultOption',3,'CancelOption',3,'Icon','warning');
        if strcmp(q,'Cancel'), return; end
        delCopies = strcmp(q,'Undo moves + delete copies');

        nOk = 0; nFail = 0;
        dlg = uiprogressdlg(fig,'Title','Undo','Value',0,'Cancelable',true);
        guard = onCleanup(@()closeIfValid(dlg));

        for i = 1:numel(m)
            if dlg.CancelRequested, log_('Undo cancelled.'); break; end
            dlg.Value   = i/numel(m);
            dlg.Message = sprintf('%d / %d', i, numel(m));
            drawnow limitrate;
            e = m(i);
            try
                if strcmp(e.op,'move')
                    if ~isfile(e.dst)
                        log_('Undo skip (not at destination): %s', e.dst); continue;
                    end
                    if isfile(e.src)
                        log_('Undo skip (source is already back): %s', e.src); continue;
                    end
                    tgtDir = fileparts(e.src);
                    if ~isfolder(tgtDir), mkdir(tgtDir); end
                    if sameVolume(e.dst, e.src)
                        [st,msg] = movefile(e.dst, e.src, 'f');
                    else
                        [st,msg] = copyfile(e.dst, e.src, 'f');
                        if st, delete(e.dst); end
                    end
                    if st
                        nOk = nOk + 1;
                    else
                        nFail = nFail + 1;
                        log_('Undo failed %s: %s', e.dst, msg);
                    end
                elseif delCopies
                    if isfile(e.dst)
                        delete(e.dst);
                        nOk = nOk + 1;
                    end
                end
            catch ME
                nFail = nFail + 1;
                log_('Undo error %s: %s', e.dst, ME.message);
            end
        end
        closeIfValid(dlg);
        log_('Undo finished: %d restored/removed, %d failed.', nOk, nFail);
    end

%% ----------------------------------------------------------------- helpers
    function log_(fmt, varargin)
        line = sprintf('%s  %s', stamp('HH:mm:ss'), sprintf(fmt, varargin{:}));
        v = txtLog.Value;
        if isscalar(v) && isempty(v{1}), v = {}; end
        v{end+1} = line;
        if numel(v) > 2000, v = v(end-1999:end); end
        txtLog.Value = v;
        scroll(txtLog,'bottom');
        drawnow limitrate;
    end

    function savePrefs
        try
            setpref('dataOrganizer','roots',S.roots);
            setpref('dataOrganizer','dst',txtDst.Value);
            setpref('dataOrganizer','patterns',txtPat.Value);
        catch
        end
    end

    function loadPrefs
        try
            if ispref('dataOrganizer','roots')
                r = getpref('dataOrganizer','roots');
                if iscell(r)
                    S.roots = r(cellfun(@(x) ischar(x) || isstring(x), r));
                    refreshRoots;
                end
            end
            if ispref('dataOrganizer','dst'),      txtDst.Value = getpref('dataOrganizer','dst'); end
            if ispref('dataOrganizer','patterns'), txtPat.Value = getpref('dataOrganizer','patterns'); end
        catch
        end
    end

    function onClose(~,~)
        savePrefs;
        delete(fig);
    end
end

%% =================================================== file-scope utilities
function s = emptyItems
s = struct('src',{},'name',{},'bytes',{},'root',{},'dst',{},'status',{},'include',{});
end

function d = skipDirs
d = {'.Trashes','.Spotlight-V100','.fseventsd','.DocumentRevisions-V100', ...
     '$RECYCLE.BIN','System Volume Information','.TemporaryItems'};
end

function tf = isHiddenPath(folder, root, name)
% True when the file name, or any folder between root and the file, is
% hidden (leading dot) or a known system folder.
tf = strncmp(name,'.',1);
if tf, return; end
rel   = relPath(fullfile(folder,'x'), root);
parts = strsplit(fileparts(rel), filesep);
sd    = skipDirs;
for i = 1:numel(parts)
    if isempty(parts{i}), continue; end
    if strncmp(parts{i},'.',1) || any(strcmp(parts{i}, sd))
        tf = true; return;
    end
end
end

function [ok,msg] = transferOne(it, isMove, allowOverwrite)
% Move or copy one file.  A move inside one volume is a rename; anything
% else is copy -> verify size -> (for a move) delete the source.
ok = false; msg = '';
tgtDir = fileparts(it.dst);
if ~isfolder(tgtDir)
    [mk,mkmsg] = mkdir(tgtDir);
    if ~mk, msg = sprintf('mkdir failed (%s)', mkmsg); return; end
end
if ~isfile(it.src),                   msg = 'source vanished';    return; end
if isfile(it.dst) && ~allowOverwrite, msg = 'destination exists'; return; end

try
    if isMove && sameVolume(it.src, it.dst)
        [st,m] = movefile(it.src, it.dst, 'f');
        if ~st, msg = m; return; end
        ok = true;
        return;
    end

    [st,m] = copyfile(it.src, it.dst, 'f');
    if ~st, msg = m; return; end
    d = dir(it.dst);
    if isempty(d) || d(1).bytes ~= it.bytes
        msg = 'size mismatch after copy - source kept';
        return;
    end
    if isMove
        delete(it.src);
        if isfile(it.src)
            msg = 'copied, but the source could not be deleted';
            return;
        end
    end
    ok = true;
catch ME
    msg = ME.message;
end
end

function n = totalBytes(items)
if isempty(items), n = 0; else, n = sum([items.bytes]); end
end

function s = stamp(fmt)
s = char(datetime('now','Format',fmt));
end

function pats = parsePatterns(val)
if isstring(val) || ischar(val), val = cellstr(val); end
txt   = strjoin(cellstr(val), ',');
parts = strtrim(strsplit(txt, {',',';',newline}));
pats  = unique(parts(~cellfun(@isempty,parts)),'stable');
end

function roots = dropNestedRoots(roots)
% Drop folders that already sit inside another listed folder.
keep = true(1,numel(roots));
for i = 1:numel(roots)
    for j = 1:numel(roots)
        if i ~= j && keep(i) && isUnder(roots{i}, roots{j})
            keep(i) = false;
        end
    end
end
roots = roots(keep);
end

function tf = isUnder(p, parent)
p      = normPath(p);
parent = normPath(parent);
if isempty(parent), tf = false; return; end
if ~endsWith(parent, filesep), parent = [parent filesep]; end
if numel(p) < numel(parent), tf = false; return; end
if ismac || ispc
    tf = strncmpi(p, parent, numel(parent));
else
    tf = strncmp(p, parent, numel(parent));
end
end

function p = normPath(p)
p = char(p);
while numel(p) > 1 && endsWith(p, filesep)
    p(end) = [];
end
end

function r = relPath(fullPath, root)
root = normPath(root);
if ~endsWith(root, filesep), root = [root filesep]; end
if numel(fullPath) > numel(root) && strncmpi(fullPath, root, numel(root))
    r = fullPath(numel(root)+1:end);
else
    [~,n,e] = fileparts(fullPath);
    r = [n e];
end
end

function b = baseName(p)
p = normPath(p);
[~,n,e] = fileparts(p);
b = [n e];
if isempty(b), b = regexprep(p,'[^\w-]','_'); end
end

function s = safeName(s)
s = strtrim(regexprep(char(s),'[<>:"/\\|?*]','_'));
if isempty(s), s = 'unnamed'; end
end

function pre = prefixOf(srcFile, depth)
parts = strsplit(fileparts(srcFile), filesep);
parts = parts(~cellfun(@isempty,parts));
if isempty(parts), pre = ''; return; end
take = parts(max(1,numel(parts)-depth+1):end);
pre  = [strjoin(cellfun(@safeName,take,'UniformOutput',false),'__') '__'];
end

function t = uniqueTarget(target, taken)
[p,n,e] = fileparts(target);
t = target;
k = 1;
while (isKey(taken, lower(t)) || isfile(t)) && k <= 9999
    t = fullfile(p, sprintf('%s_%d%s', n, k, e));
    k = k + 1;
end
end

function tf = isLink(p)
% True when the last path component itself is a symlink / alias.
tf = false;
try
    f      = java.io.File(p);
    parent = char(f.getParent());
    if isempty(parent), return; end
    expect = fullfile(char(java.io.File(parent).getCanonicalPath), char(f.getName()));
    tf     = ~strcmp(expect, char(f.getCanonicalPath));
catch
end
end

function tf = sameVolume(a, b)
% macOS / Windows aware guess at whether two paths live on one volume.
tf = strcmpi(volumeOf(a), volumeOf(b));
end

function v = volumeOf(p)
p = char(p);
if ispc
    v = upper(p(1:min(2,numel(p))));
    return;
end
tok = regexp(p,'^/Volumes/[^/]+','match','once');
if isempty(tok), v = '/'; else, v = tok; end
end

function s = humanBytes(b)
if isempty(b) || ~isfinite(b), s = '-'; return; end
u = {'B','KB','MB','GB','TB','PB'};
k = 1;
while b >= 1024 && k < numel(u)
    b = b/1024;
    k = k + 1;
end
if k == 1, s = sprintf('%d %s', round(b), u{k});
else,      s = sprintf('%.1f %s', b, u{k});
end
end

function s = humanTime(sec)
if sec < 60
    s = sprintf('%.1f s', sec);
elseif sec < 3600
    s = sprintf('%d m %02d s', floor(sec/60), round(mod(sec,60)));
else
    s = sprintf('%d h %02d m', floor(sec/3600), round(mod(sec,3600)/60));
end
end

function d = uigetdir_(startDir, title)
if isempty(startDir) || ~isfolder(startDir), startDir = pwd; end
d = uigetdir(startDir, title);
if isequal(d,0), d = ''; end
end

function closeIfValid(dlg)
try
    if isvalid(dlg), close(dlg); end
catch
end
end

%% ===================================================== self-test (no GUI)
function selfTest
fprintf('--- dataOrganizer self-test ---\n');
nFail = 0;

% ---- pure helpers
check(isequal(parsePatterns({'frame1.bin, frame2.bin','output_data.mat'}), ...
              {'frame1.bin','frame2.bin','output_data.mat'}), 'parsePatterns splits commas and lines');
check(isequal(parsePatterns({'a.bin,,  ,a.bin'}), {'a.bin'}), 'parsePatterns drops blanks and duplicates');
check(strcmp(humanBytes(0),'0 B') && strcmp(humanBytes(1536),'1.5 KB') && ...
      strcmp(humanBytes(3*1024^3),'3.0 GB'), 'humanBytes formats sizes');
check(isUnder('/a/b/c.bin','/a/b') && ~isUnder('/a/bb/c.bin','/a/b') && ...
      ~isUnder('/a/b','/a/b'), 'isUnder handles prefixes without false positives');
check(isequal(dropNestedRoots({'/a','/a/b','/c'}), {'/a','/c'}), 'dropNestedRoots removes nested folders');
check(strcmp(relPath(fullfile('/r','s1','s2','f.bin'), '/r'), ...
             fullfile('s1','s2','f.bin')), 'relPath is relative to the root');
check(strcmp(relPath('/other/f.bin','/r'),'f.bin'), 'relPath falls back to the bare name');
check(strcmp(prefixOf(fullfile('/r','expA','run3','f.bin'),1),'run3__'), 'prefixOf depth 1');
check(strcmp(prefixOf(fullfile('/r','expA','run3','f.bin'),2),'expA__run3__'), 'prefixOf depth 2');
check(strcmp(safeName('a/b:c'),'a_b_c'), 'safeName strips path characters');
check(sameVolume('/Volumes/D1/a','/Volumes/D1/b') && ...
      ~sameVolume('/Volumes/D1/a','/Volumes/D2/b') && ...
      ~sameVolume('/Users/x/a','/Volumes/D1/b'), 'sameVolume compares macOS volumes');

emptyMap = containers.Map('KeyType','char','ValueType','logical');
check(strcmp(uniqueTarget('/nope/f.bin', emptyMap), '/nope/f.bin'), 'uniqueTarget keeps a free name');
usedMap = containers.Map({'/nope/f.bin'}, {true});
check(strcmp(uniqueTarget('/nope/f.bin', usedMap), fullfile('/nope','f_1.bin')), ...
      'uniqueTarget renames a taken name');

% ---- real files in a temporary tree
tmp  = tempname;
src  = fullfile(tmp,'sourceDisk');
dst  = fullfile(tmp,'collected');
mkdir(fullfile(src,'expA','run1'));
mkdir(fullfile(src,'expA','run2'));
mkdir(fullfile(src,'.hidden','run3'));
cleanup = onCleanup(@() rmdirIfPresent(tmp));

writeBytes(fullfile(src,'expA','run1','frame1.bin'), 1000);
writeBytes(fullfile(src,'expA','run1','frame2.bin'), 2000);
writeBytes(fullfile(src,'expA','run1','output_data.mat'), 300);
writeBytes(fullfile(src,'expA','run2','frame1.bin'), 1100);
writeBytes(fullfile(src,'expA','run2','ignore_me.bin'), 10);
writeBytes(fullfile(src,'.hidden','run3','frame1.bin'), 7);

pats  = {'frame1.bin','frame2.bin','output_data.mat'};
found = {};
for ip = 1:numel(pats)
    d = dir(fullfile(src,'**',pats{ip}));
    d = d(~[d.isdir]);
    for k = 1:numel(d)
        if isHiddenPath(d(k).folder, src, d(k).name), continue; end
        found{end+1} = fullfile(d(k).folder, d(k).name); %#ok<AGROW>
    end
end
check(numel(found) == 4, sprintf('recursive scan found the 4 wanted files (got %d)', numel(found)));
check(~any(contains(found,'ignore_me')), 'scan ignores unlisted names');
check(~any(contains(found,[filesep '.hidden' filesep])), 'scan skips hidden folders');

% mirror layout, two files with the same name in different subfolders
taken = containers.Map('KeyType','char','ValueType','logical');
nMoved = 0;
for i = 1:numel(found)
    rel    = relPath(found{i}, src);
    target = fullfile(dst, safeName(baseName(src)), rel);
    target = uniqueTarget(target, taken);
    taken(lower(target)) = true;
    it = struct('src',found{i},'dst',target,'bytes',dir(found{i}).bytes);
    [ok,msg] = transferOne(it, true, false);
    check(ok, sprintf('moved %s%s', rel, ternary(ok,'',[' (' msg ')'])));
    if ok, nMoved = nMoved + 1; end
end
check(nMoved == 4, 'all four files were moved');
check(isfile(fullfile(dst,'sourceDisk','expA','run1','frame1.bin')) && ...
      isfile(fullfile(dst,'sourceDisk','expA','run2','frame1.bin')), ...
      'mirror layout keeps same-named files apart');
check(dir(fullfile(dst,'sourceDisk','expA','run1','frame2.bin')).bytes == 2000, ...
      'moved file keeps its size');
check(~isfile(fullfile(src,'expA','run1','frame1.bin')), 'source is gone after a move');
check(isfile(fullfile(src,'expA','run2','ignore_me.bin')), 'unlisted file was left alone');

% collision guard: a second move onto an existing destination must refuse
writeBytes(fullfile(src,'expA','run1','frame1.bin'), 999);
it = struct('src',fullfile(src,'expA','run1','frame1.bin'), ...
            'dst',fullfile(dst,'sourceDisk','expA','run1','frame1.bin'),'bytes',999);
[ok,msg] = transferOne(it, true, false);
check(~ok && strcmp(msg,'destination exists'), 'refuses to overwrite unless allowed');
check(dir(it.dst).bytes == 1000, 'the existing destination file is untouched');
[ok,~] = transferOne(it, true, true);
check(ok && dir(it.dst).bytes == 999, 'overwrites when explicitly allowed');

% flat + prefix layout
writeBytes(fullfile(src,'expA','run1','frame1.bin'), 1000);
writeBytes(fullfile(src,'expA','run2','frame1.bin'), 1100);
flat  = fullfile(tmp,'flat');
taken = containers.Map('KeyType','char','ValueType','logical');
for p = {fullfile(src,'expA','run1','frame1.bin'), fullfile(src,'expA','run2','frame1.bin')}
    target = fullfile(flat, [prefixOf(p{1},1) 'frame1.bin']);
    target = uniqueTarget(target, taken);
    taken(lower(target)) = true;
    transferOne(struct('src',p{1},'dst',target,'bytes',dir(p{1}).bytes), false, false);
end
check(isfile(fullfile(flat,'run1__frame1.bin')) && isfile(fullfile(flat,'run2__frame1.bin')), ...
      'flat + prefix layout produces unique names');
check(isfile(fullfile(src,'expA','run1','frame1.bin')), 'copy leaves the source in place');

% manifest round trip
mani = struct('op',{'move','move'},'src',{'/a/f1.bin','/a/f2.bin'}, ...
              'dst',{'/b/f1.bin','/b/f2.bin'},'bytes',{10,20});
back = jsondecode(jsonencode(mani,'PrettyPrint',true));
check(numel(back) == 2 && strcmp(back(2).src,'/a/f2.bin'), 'manifest survives a JSON round trip');

fprintf('--- %s ---\n', ternary(nFail == 0, 'all checks passed', sprintf('%d CHECK(S) FAILED', nFail)));
if nFail > 0, error('dataOrganizer:selfTestFailed','%d check(s) failed', nFail); end

    function check(cond, what)
        if cond
            fprintf('  ok    %s\n', what);
        else
            fprintf('  FAIL  %s\n', what);
            nFail = nFail + 1;
        end
    end
end

function writeBytes(p, n)
fid = fopen(p,'w');
fwrite(fid, zeros(n,1,'uint8'));
fclose(fid);
end

function rmdirIfPresent(p)
if isfolder(p), rmdir(p,'s'); end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
