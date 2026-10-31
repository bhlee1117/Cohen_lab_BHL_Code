function mc_log(cfg, level, fmt, varargin)
%MC_LOG  Timestamped progress logging with optional phone notification.
%
%   mc_log(cfg, level, fmt, ...) prints a timestamped message to the command
%   window, appends it to cfg.logfile, and - for level >= 1 - posts it to
%   cfg.notify_url so it shows up on your phone.
%
%   cfg fields
%     .logfile       full path of the text log (created if missing)
%     .notify_url    webhook to POST to; '' disables notifications
%                      ntfy  : 'https://ntfy.sh/<your-private-topic>'
%                      Slack : 'https://hooks.slack.com/services/...'
%     .notify_title  short title shown on the notification (ntfy only)
%
%   level
%     0  command window + log file only
%     1  quiet notification (ntfy priority "min": shows up, no sound)
%     2  normal notification
%     3  urgent notification (ntfy priority "high")
%
%   Notifications are best effort: a failed POST is reported but never
%   interrupts the analysis.
%
%   Example
%     cfg.logfile     = '~/Documents/MC_progress.log';
%     cfg.notify_url  = 'https://ntfy.sh/bhl-matlab-7f3a2c91';
%     cfg.notify_title= 'MATLAB motion correction';
%     mc_log(cfg, 2, 'DONE file#%d', 20);

msg  = sprintf(fmt, varargin{:});
line = sprintf('%s  %s', char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), msg);

fprintf('%s\n', line);

% ---- append to the log file -----------------------------------------
if isfield(cfg,'logfile') && ~isempty(cfg.logfile)
    fid = fopen(cfg.logfile, 'a');
    if fid > 0
        fprintf(fid, '%s\n', line);
        fclose(fid);
    else
        warning('mc_log:noLogFile','Could not open log file %s', cfg.logfile);
    end
end

% ---- push notification ----------------------------------------------
if level < 1 || ~isfield(cfg,'notify_url') || isempty(cfg.notify_url)
    return
end

ntitle = 'MATLAB';
if isfield(cfg,'notify_title') && ~isempty(cfg.notify_title)
    ntitle = cfg.notify_title;
end

try
    if contains(cfg.notify_url, 'hooks.slack.com')
        % Slack incoming webhook expects a JSON payload
        webwrite(cfg.notify_url, struct('text', sprintf('%s: %s', ntitle, msg)), ...
            weboptions('MediaType','application/json','Timeout',15));
    else
        % ntfy.sh and compatible: plain-text body, metadata in headers
        prio = {'min','default','high'};
        hdrs = {'Title', ntitle; 'Priority', prio{min(level,3)}};
        webwrite(cfg.notify_url, msg, ...
            weboptions('MediaType','text/plain','Timeout',15,'HeaderFields',hdrs));
    end
catch ME
    fprintf('  (notification failed: %s)\n', ME.message);
end
end
