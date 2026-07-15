function windowsPath = macToWindowsPath(macPath, server)
% MACTOWINDOWSPATH Convert a Mac-style /Volumes/<share>/... path to a
% Windows UNC path \\server\<share>\...
% If the input does not start with /Volumes/, it is returned unchanged.
%
%   windowsPath = macToWindowsPath(macPath)
%   windowsPath = macToWindowsPath(macPath, server)
%
%   Example:
%       p = "/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Prism_V2+Glu_Data_Arrangement.xlsx";
%       macToWindowsPath(p)
%
%   ans =
%       \\b-nfs01-p.rc.fas.harvard.edu\cohen_lab\Lab\Labmembers\Byung Hun Lee\Data\Prism_V2+Glu_Data_Arrangement.xlsx

    if nargin < 2
        server = "b-nfs01-p.rc.fas.harvard.edu";
    end

    macPath = string(macPath);

    % Match and strip the leading /Volumes/ part
    tok = regexp(macPath, "^/Volumes/(.+)$", "tokens", "once");

    if isempty(tok)
        % Not a /Volumes/ path (e.g. already a UNC path, or a different
        % format) -- leave it unchanged
        windowsPath = macPath;
        return
    end

    remainder = tok{1};

    % Convert forward slashes to backslashes
    remainderWin = strrep(remainder, "/", "\");

    % Build the UNC path
    windowsPath = "\\" + server + "\" + remainderWin;
end