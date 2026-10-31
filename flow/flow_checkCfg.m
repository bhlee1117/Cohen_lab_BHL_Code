function [cfg, rep] = flow_checkCfg(cfg, quiet)
% [cfg, rep] = flow_checkCfg(cfg)
% [cfg, rep] = flow_checkCfg(cfg, true)      % fix it, say nothing
%
% Reconcile a config with the CURRENT flow_defaults: fill in knobs that have been added
% since it was made, drop knobs that have been removed, and say what changed.
%
% WHY THIS EXISTS. data goes through flow_checkData; cfg went through nothing. A config
% saved before a knob existed passes straight into the pipeline and fails at the first
% line that reads the missing field -- deep inside a stage, as "Unrecognized field name",
% with nothing to say that the real problem is an old cfg rather than a bug. That has
% happened for every knob added since the package was written. This turns the whole class
% of failure into one message at the front door.
%
% WHAT IT DOES NOT DO: guess intent across a RENAME. If voltLo/voltHi/voltMin became
% voltCut, this adds voltCut at its default and reports the three it dropped -- it does
% not try to translate the old values into the new one, because the mapping is a judgement
% about what the knob now means and only the author of the change knows it. Read the
% dropped list and re-set anything you had deliberately tuned.
%
% THE STAMPS WILL MOVE. Adding or removing a field changes what flow_stamp hashes, so a
% cfg that needed fixing will mark its stage -- and everything downstream -- STALE against
% results computed before. That is honest: those results were computed with a different
% parameter set. Field ORDER does not matter (flow_stamp sorts), only the set of names.
%
% INPUT   cfg    a config struct; [] or omitted gives flow_defaults()
%         quiet  true to suppress the report (default false)
% OUTPUT  cfg    the reconciled config
%         rep    .added   {'stage.field', ...} filled in from defaults
%                .removed {'stage.field', ...} dropped as unknown
%                .stages  stages copied in whole because they were missing
%                .ok      true when nothing had to change
%
% 2026.08.30, Byung Hun Lee

if nargin < 2 || isempty(quiet), quiet = false; end
if nargin < 1 || isempty(cfg),   cfg = flow_defaults(); end
assert(isstruct(cfg), 'cfg must be a struct.');

dflt = flow_defaults();
rep = struct('added', {{}}, 'removed', {{}}, 'stages', {{}}, 'ok', true);

for s = fieldnames(dflt)'
    st = s{1};
    if ~isfield(cfg, st) || ~isstruct(cfg.(st))
        cfg.(st) = dflt.(st);
        rep.stages{end+1} = st;                                     %#ok<AGROW>
        continue
    end
    for g = setdiff(fieldnames(dflt.(st)), fieldnames(cfg.(st)))'
        cfg.(st).(g{1}) = dflt.(st).(g{1});
        rep.added{end+1} = [st '.' g{1}];                           %#ok<AGROW>
    end
    for g = setdiff(fieldnames(cfg.(st)), fieldnames(dflt.(st)))'
        cfg.(st) = rmfield(cfg.(st), g{1});
        rep.removed{end+1} = [st '.' g{1}];                         %#ok<AGROW>
    end
end
% whole stages that no longer exist
for s = setdiff(fieldnames(cfg), fieldnames(dflt))'
    cfg = rmfield(cfg, s{1});
    rep.removed{end+1} = s{1};                                      %#ok<AGROW>
end

rep.ok = isempty(rep.added) && isempty(rep.removed) && isempty(rep.stages);
if quiet || rep.ok, return, end

fprintf('[flow] cfg reconciled with the current flow_defaults:\n');
if ~isempty(rep.stages)
    fprintf('       %d stage(s) missing entirely, copied from defaults: %s\n', ...
        numel(rep.stages), strjoin(rep.stages, ', '));
end
if ~isempty(rep.added)
    fprintf('       %d knob(s) ADDED at their default value:\n         %s\n', ...
        numel(rep.added), strjoin(rep.added, ', '));
end
if ~isempty(rep.removed)
    fprintf('       %d knob(s) REMOVED as unknown to this version:\n         %s\n', ...
        numel(rep.removed), strjoin(rep.removed, ', '));
    fprintf(['       If one of those was RENAMED rather than deleted, its replacement is ' ...
             'now sitting at\n       the default -- re-set it if you had tuned the old ' ...
             'one.\n']);
end
fprintf(['       Stages whose parameter SET changed will read STALE against older ' ...
         'results, because\n       they were computed with a different set.\n']);
end
