function info = flow_paramInfo()
% info = flow_paramInfo()
%
% What every stage and every parameter means, for the GUI's note column.
%
% The parameter texts are PARSED OUT OF flow_defaults.m rather than written again here.
% Duplicating them would guarantee the two drift apart -- the comment next to the value
% is the one a person edits when they change what a knob does, so that comment is the
% single source of truth and this function just reads it.
%
% OUTPUT  info
%   .stage.(name)  one line describing what that stage does
%   .param.(stage).(name)  the comment from flow_defaults, flattened to one line
%
% If flow_defaults is reformatted so the parse fails, descriptions come back empty and
% the GUI simply shows no note -- it never guesses.
%
% 2026.08.28, Byung Hun Lee

info.stage = struct( ...
 'bandpass',   ['Spatial band-pass of the motion-residual movie. Keeps only structure ' ...
                'on the vessel''s spatial scale (~6-30 um), which isolates blood using ' ...
                'no temporal assumption at all. Optionally projects the measured ' ...
                'voltage out first. -> S'], ...
 'track',      ['Detect red blood cells in S frame by frame and link them into ' ...
                'tracklets. Gives the flow SPEED and direction per vessel, which is ' ...
                'what sets the frequency band later (f = speed / spacing). -> tracklets'], ...
 'segments',   ['Turn the tracklets into VESSEL SEGMENTS: a density map of where cells ' ...
                'were actually found, thresholded, split into connected patches, and ' ...
                'each patch matched to the tracklets inside it. Every later stage is ' ...
                'built on these masks. -> segment masks'], ...
 'regressors', ['Per segment, build a phase-continuous oscillator from the measured ' ...
                'signal (band-pass + Hilbert, or complex PCA) and assemble the design ' ...
                'matrix. Scores each row and drops the ones that look like noise or ' ...
                'like the neuron. -> design matrix'], ...
 'apply',      ['Fit that design to the whole recording and subtract it. Re-derives the ' ...
                'modes on the full record, checks whether the band drifted, and tests ' ...
                'whether a time-varying coefficient does better. -> corrected movie'], ...
 'compare',    ['What the correction does to the extracted TRACES: how much variance it ' ...
                'removes and where, in which frequencies, and whether the spikes ' ...
                'survive it.']);

info.param = struct();
fp = which('flow_defaults');
if isempty(fp), return, end
txt = strsplit(fileread(fp), newline);

% Parameters are collected into `cur` as they are read, and moved into info.param once
% the section's closing 'cfg.<stage> = c;' line says which stage they belonged to.
cur = struct();  pend = '';
for k = 1:numel(txt)
    L = txt{k};
    % 'cfg.<stage> = c;' closes a section and names it
    tok = regexp(L, '^\s*cfg\.(\w+)\s*=\s*c\s*;', 'tokens', 'once');
    if ~isempty(tok)
        info.param.(tok{1}) = cur;
        cur = struct();  pend = '';  continue
    end
    % 'c.<name> = <value>;   % description'   (the value may itself contain a % in a string)
    tok = regexp(L, '^\s*c\.(\w+)\s*=[^%]*%\s*(.*)$', 'tokens', 'once');
    if ~isempty(tok)
        pend = tok{1};  cur.(pend) = strtrim(tok{2});  continue
    end
    % a value line with no trailing comment still opens a parameter
    tok = regexp(L, '^\s*c\.(\w+)\s*=', 'tokens', 'once');
    if ~isempty(tok)
        pend = tok{1};  cur.(pend) = '';  continue
    end
    % a comment-only continuation line extends the parameter opened above
    if ~isempty(pend)
        tok = regexp(L, '^\s*%\s?(.*)$', 'tokens', 'once');
        if ~isempty(tok) && ~isempty(strtrim(tok{1})) && ...
                isempty(regexp(tok{1}, '^[=-]{3,}', 'once'))
            cur.(pend) = strtrim([cur.(pend) ' ' tok{1}]);
            continue
        end
    end
    pend = '';                       % anything else ends the continuation
end
end
