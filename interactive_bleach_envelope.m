function [env, info, hp] = interactive_bleach_envelope(trace, fs, thr0, hpWin0, binOpts)
% INTERACTIVE_BLEACH_ENVELOPE  Spike-height envelope (photobleaching) of a voltage trace.
%   [env, info] = interactive_bleach_envelope(trace, fs, thr0)            % trace already high-passed
%   [env, info, hp] = interactive_bleach_envelope(trace, fs, thr0, hpWin0) % high-pass set here
%   [env, info, hp] = interactive_bleach_envelope(trace, fs, thr0, hpWin0, binOpts)
%
%   With hpWin0, trace is the un-filtered (noise-normalised) trace and the window sets the
%   high-pass first: hp = trace - movmedian(trace, hpWin) (frames). Changing it re-runs the
%   spike finding and the envelope fit.
%
%   Putative spikes are found with find_spike_bh(trace, thr1, thr2) (thresholds set in the
%   window), their heights are binned in time and each bin is summarised by a percentile of the
%   heights ('Percentile', default 75; bins with fewer than 'Min spikes/bin' spikes, default 1,
%   are skipped; at least 4 bins are needed, otherwise no correction) and an envelope is fitted:
%       'exp'    : a*exp(-t/tau) + c
%       'exp2'   : a1*exp(-t/tau1) + a2*exp(-t/tau2) + c
%       'smooth' : moving average of the bin values, interpolated
%   env is the envelope normalised to 1 at the start of the trace; trace ./ env has a flat
%   spike height in the units of the start of the recording.
%
%   INPUT
%     trace  : 1 x T voltage trace (spikes positive-going); high-passed already unless hpWin0
%     fs     : frame rate (Hz), default 1 (time in frames)
%     thr0   : initial [thr1 thr2], default [5 2]
%     hpWin0 : initial high-pass window (frames); [] = trace is already high-passed
%     binOpts: optional struct with initial .minSpikes (default 1) and .pct (default 75)
%   OUTPUT
%     env   : 1 x T normalised envelope (ones when 'No correction' is ticked)
%     info  : .thr, .method, .bin_s, .smooth_s, .minSpikes, .pct, .params, .nSpikes, .corrected, .hpWin
%     hp    : 1 x T high-passed trace the envelope was fitted on

    if nargin < 2 || isempty(fs), fs = 1; end
    if nargin < 3 || isempty(thr0), thr0 = [5 2]; end
    if nargin < 4, hpWin0 = []; end
    minSpikes0 = 1; pct0 = 75;
    if nargin > 4 && isstruct(binOpts)
        if isfield(binOpts, 'minSpikes') && ~isempty(binOpts.minSpikes), minSpikes0 = binOpts.minSpikes; end
        if isfield(binOpts, 'pct') && ~isempty(binOpts.pct), pct0 = binOpts.pct; end
    end
    raw = double(trace(:)'); raw(~isfinite(raw)) = 0;
    useHP = ~isempty(hpWin0);
    hpWin = NaN; trace = raw; hpDone = NaN;
    T = numel(trace); t = (0:T-1) / fs;
    if fs == 1, tUnit = 'frame'; else, tUnit = 's'; end
    methods = {'exp','exp2','smooth'};

    env = ones(1,T); info = struct();

    f = figure('Name','Interactive Bleach Envelope','NumberTitle','off','Position',[100 100 1400 800]);
    try, theme(f, 'light'); catch, end

    ax1 = axes('Parent',f,'Position',[0.05 0.62 0.92 0.33]);
    hTr = plot(ax1, t, trace, 'k'); hold(ax1,'on'); grid(ax1,'on');
    hSp = plot(ax1, nan, nan, 'ro', 'MarkerSize', 5);
    hEnv1 = plot(ax1, nan, nan, 'b-', 'LineWidth', 2);
    ylabel(ax1,'Amplitude'); xlabel(ax1, tUnit);
    legend(ax1, {'Trace','Putative spikes','Envelope'}, 'Location','northeast');

    ax2 = axes('Parent',f,'Position',[0.05 0.25 0.44 0.28]);
    hH  = plot(ax2, nan, nan, '.', 'Color',[0.6 0.6 0.6]); hold(ax2,'on'); grid(ax2,'on');
    hB  = plot(ax2, nan, nan, 'bs', 'MarkerFaceColor','b', 'MarkerSize',5);
    hEnv2 = plot(ax2, nan, nan, 'r-', 'LineWidth', 2);
    ylabel(ax2,'spike height'); xlabel(ax2, tUnit); title(ax2,'putative spike heights, bin percentile (blue), envelope','FontSize',9);

    ax3 = axes('Parent',f,'Position',[0.55 0.25 0.42 0.28]);
    hHc = plot(ax3, nan, nan, '.', 'Color',[0.6 0.6 0.6]); hold(ax3,'on'); grid(ax3,'on');
    hBc = plot(ax3, nan, nan, 'bs', 'MarkerFaceColor','b', 'MarkerSize',5);
    ylabel(ax3,'corrected height'); xlabel(ax3, tUnit); title(ax3,'after dividing by the envelope (should be flat)','FontSize',9);

    % controls
    y = 0.1; h = 0.04;
    uicontrol('Style','text','String','High-pass window (frames):','Units','normalized','Position',[0.03 y+0.06 0.12 h],'HorizontalAlignment','right');
    hpEdit = uicontrol('Style','edit','String',num2str(hpWin0),'Units','normalized','Position',[0.16 y+0.06 0.05 h],'Callback',@update);
    hpMs = uicontrol('Style','text','String','','Units','normalized','Position',[0.22 y+0.06 0.3 h],'HorizontalAlignment','left');
    if ~useHP, set(hpEdit,'Enable','off','String','-'); set(hpMs,'String','trace given already high-passed'); end
    uicontrol('Style','text','String','Min spikes/bin:','Units','normalized','Position',[0.54 y+0.06 0.08 h],'HorizontalAlignment','right');
    minEdit = uicontrol('Style','edit','String',num2str(minSpikes0),'Units','normalized','Position',[0.63 y+0.06 0.05 h],'Callback',@update);
    uicontrol('Style','text','String','Percentile:','Units','normalized','Position',[0.69 y+0.06 0.07 h],'HorizontalAlignment','right');
    pctEdit = uicontrol('Style','edit','String',num2str(pct0),'Units','normalized','Position',[0.77 y+0.06 0.05 h],'Callback',@update);
    uicontrol('Style','text','String','Threshold1:','Units','normalized','Position',[0.03 y 0.06 h],'HorizontalAlignment','right');
    th1Edit = uicontrol('Style','edit','String',num2str(thr0(1)),'Units','normalized','Position',[0.10 y 0.05 h],'Callback',@update);
    uicontrol('Style','text','String','Threshold2:','Units','normalized','Position',[0.16 y 0.06 h],'HorizontalAlignment','right');
    th2Edit = uicontrol('Style','edit','String',num2str(thr0(2)),'Units','normalized','Position',[0.23 y 0.05 h],'Callback',@update);
    uicontrol('Style','text','String',sprintf('Bin (%s):',tUnit),'Units','normalized','Position',[0.29 y 0.06 h],'HorizontalAlignment','right');
    binEdit = uicontrol('Style','edit','String',num2str(round(20*fs)/fs),'Units','normalized','Position',[0.36 y 0.05 h],'Callback',@update);
    uicontrol('Style','text','String','Method:','Units','normalized','Position',[0.42 y 0.05 h],'HorizontalAlignment','right');
    methPop = uicontrol('Style','popupmenu','String',methods,'Value',1,'Units','normalized','Position',[0.48 y 0.07 h],'Callback',@update);
    uicontrol('Style','text','String',sprintf('Smooth (%s):',tUnit),'Units','normalized','Position',[0.56 y 0.07 h],'HorizontalAlignment','right');
    smEdit = uicontrol('Style','edit','String',num2str(round(60*fs)/fs),'Units','normalized','Position',[0.64 y 0.05 h],'Callback',@update);
    noCorr = uicontrol('Style','checkbox','String','No correction','Units','normalized','Position',[0.70 y 0.08 h],'Callback',@update);
    uicontrol('Style','pushbutton','String','Update','Units','normalized','Position',[0.79 y 0.08 h],'Callback',@update);
    uicontrol('Style','pushbutton','String','Done','Units','normalized','Position',[0.89 y 0.08 h],'Callback','uiresume(gcbf);');
    infoTxt = uicontrol('Style','text','String','','Units','normalized','Position',[0.03 0.03 0.9 h],'HorizontalAlignment','left');

    update();
    uiwait(f);
    if isvalid(f), update(); close(f); end

    % ================= nested =================
    function update(~, ~)
        if useHP
            hpWin = max(3, round(str2double(hpEdit.String)));
            if ~isequal(hpWin, hpDone)
                trace = raw - movmedian(raw, hpWin, 2);
                set(hTr, 'YData', trace); hpDone = hpWin;
            end
            set(hpMs, 'String', sprintf('= %.0f ms  (trace - movmedian(trace, window))', 1000*hpWin/fs));
        end
        thr = [str2double(th1Edit.String) str2double(th2Edit.String)];
        bin = str2double(binEdit.String) * fs;                     % frames
        minN = max(1, round(str2double(minEdit.String)));          % spikes needed in a bin
        pct = min(100, max(0, str2double(pctEdit.String)));        % bin value: this percentile of heights
        sm  = str2double(smEdit.String) * fs;
        meth = methods{methPop.Value};

        loc = find(find_spike_bh(trace, thr(1), thr(2)));
        ht  = trace(loc);
        set(hSp, 'XData', t(loc), 'YData', ht);
        set(hH,  'XData', t(loc), 'YData', ht);

        % bin values: percentile pct of the heights (bins with fewer than minN spikes are not used)
        edges = 0:bin:T+bin; nb = numel(edges)-1;
        [tb, hb, nbc] = deal(nan(1,nb));
        for k = 1:nb
            in = loc > edges(k) & loc <= edges(k+1);
            nbc(k) = sum(in);
            if nbc(k) >= minN, tb(k) = median(loc(in)); hb(k) = prctile(ht(in), pct); end
        end
        ok = ~isnan(hb); tb = tb(ok); hb = hb(ok); w = sqrt(nbc(ok));
        set(hB, 'XData', (tb-1)/fs, 'YData', hb);

        params = []; model = ones(1,T);
        if numel(hb) >= 4 && ~noCorr.Value
            x = 1:T;
            switch meth
                case 'exp'
                    p0 = [log(max(hb(1)-hb(end),eps)), log(T/3), hb(end)];
                    fun = @(p,xx) exp(p(1))*exp(-xx/exp(p(2))) + p(3);
                    params = fminsearch(@(p) sum(w.*(fun(p,tb)-hb).^2), p0, optimset('Display','off','MaxFunEvals',4000,'MaxIter',4000));
                    model = fun(params, x);
                    params = struct('a',exp(params(1)), 'tau_s',exp(params(2))/fs, 'c',params(3));
                case 'exp2'
                    d = max(hb(1)-hb(end), eps);
                    p0 = [log(d/2), log(T/20), log(d/2), log(T/2), hb(end)];
                    fun = @(p,xx) exp(p(1))*exp(-xx/exp(p(2))) + exp(p(3))*exp(-xx/exp(p(4))) + p(5);
                    params = fminsearch(@(p) sum(w.*(fun(p,tb)-hb).^2), p0, optimset('Display','off','MaxFunEvals',1e4,'MaxIter',1e4));
                    model = fun(params, x);
                    params = struct('a1',exp(params(1)), 'tau1_s',exp(params(2))/fs, 'a2',exp(params(3)), ...
                        'tau2_s',exp(params(4))/fs, 'c',params(5));
                case 'smooth'
                    nsm = max(1, round(sm / bin));
                    hs = movmean(hb, nsm);
                    model = interp1(tb, hs, x, 'linear', 'extrap');
                    model(x < tb(1)) = hs(1); model(x > tb(end)) = hs(end);
                    params = struct('tb_s',(tb-1)/fs, 'hb_smooth',hs);
            end
            model = max(model, 0.05*max(model));                  % guard against blow-up
            env = model / model(1);
        else
            env = ones(1,T);
        end
        set(hEnv2, 'XData', t, 'YData', env * model(1));
        set(hEnv1, 'XData', t, 'YData', env * model(1));

        % corrected heights
        hc = ht ./ env(loc);
        set(hHc, 'XData', t(loc), 'YData', hc);
        hbc = nan(1,nb); tbc = nan(1,nb);
        for k = 1:nb
            in = loc > edges(k) & loc <= edges(k+1);
            if sum(in) >= minN, tbc(k) = median(loc(in)); hbc(k) = prctile(hc(in), pct); end
        end
        set(hBc, 'XData', (tbc-1)/fs, 'YData', hbc);

        info = struct('thr',thr, 'method',meth, 'bin_s',bin/fs, 'smooth_s',sm/fs, 'minSpikes',minN, 'pct',pct, 'params',params, ...
            'nSpikes',numel(loc), 'corrected',~noCorr.Value && numel(hb) >= 4, 'hpWin',hpWin);
        hp = trace;
        set(infoTxt, 'String', sprintf('%d putative spikes, %d bins with >= %d spikes | envelope end/start = %.2f | corrected bin P%g: first %.2f, last %.2f', ...
            numel(loc), numel(hb), minN, env(end), pct, hbc(find(~isnan(hbc),1)), hbc(find(~isnan(hbc),1,'last'))));
        drawnow;
    end
end
