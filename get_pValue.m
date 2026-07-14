function p_values = get_pValue(data, isPaired, varargin)
% get_pValue - pairwise p-values for paired or unpaired groups with selectable test.
%
% p_values = get_pValue(data, isPaired)
% p_values = get_pValue(..., 'Method','ttest2')      % unpaired
% p_values = get_pValue(..., 'Method','kstest2')     % unpaired
% p_values = get_pValue(..., 'Method','ttest')       % paired
% p_values = get_pValue(..., 'Method','signrank')    % paired
% p_values = get_pValue(..., 'TestFcn',@(x,y) ...)   % custom function -> returns p
%
% data:   numeric matrix (paired) or cell array (unpaired)
% isPaired: logical
%
% Options:
%   'Method'  : char/string, default 'auto'
%   'TestFcn' : function handle (x,y)->p. Overrides Method if provided.

% Allow method to be passed directly: get_pValue(data, isPaired, 'ttest')
if numel(varargin) == 1 && (ischar(varargin{1}) || isstring(varargin{1}))
    varargin = {'Method', varargin{1}};
end

ip = inputParser;
ip.addParameter('Method','auto',@(s)ischar(s)||isstring(s));
ip.addParameter('TestFcn',[],@(f)isempty(f)||isa(f,'function_handle'));
ip.parse(varargin{:});
method  = lower(string(ip.Results.Method));
testFcn = ip.Results.TestFcn;

if isPaired
    numGroups = size(data, 2);
else
    numGroups = numel(data);
end

p_values = nan(numGroups);
printedHeader = false;

for i = 1:numGroups
    for j = i+1:numGroups

        % ----- pull clean vectors -----
        if isPaired
            pairedData = [data(:,i), data(:,j)];
            mask = all(~isnan(pairedData), 2);
            x = pairedData(mask,1);
            y = pairedData(mask,2);
        else
            x = data{i}; x = x(~isnan(x));
            y = data{j}; y = y(~isnan(y));
        end

        % Skip comparisons where either group has no data (e.g. an empty bin);
        % ranksum/ttest error on empty input. Report NaN for that pair.
        if isempty(x) || isempty(y)
            p_values(i,j) = NaN; p_values(j,i) = NaN;
            continue
        end

        % ----- run requested test -----
        if ~isempty(testFcn)
            pval = testFcn(x,y);

        else
            if method == "auto"
                if isPaired
                    method = "ttest";
                else
                    method = "ranksum";
                end
            end

            switch method
                % Paired
                case "ttest"
                    [~,pval] = ttest(x,y);
                case "signrank"
                    pval = signrank(x,y);
                case "kstest"   % paired KS on differences (less common)
                    pval = kstest(x - y);

                % Unpaired
                case "ttest2"
                    [~,pval] = ttest2(x,y);
                case "ranksum"
                    pval = ranksum(x,y);
                case "kstest2"
                    [~,pval] = kstest2(x,y);

                otherwise
                    error("Unknown Method '%s'. Try: ttest, signrank, ttest2, ranksum, kstest2, or supply 'TestFcn'.", method);
            end
        end

        p_values(i,j) = pval;
        p_values(j,i) = pval;

        % Print header once after method is resolved
        if ~printedHeader
            if ~isempty(testFcn)
                testLabel = 'custom TestFcn';
            else
                testLabel = char(method);
            end
            fprintf('\n--- Statistical Tests: %s ---\n', testLabel);
            printedHeader = true;
        end
        fprintf('  Group %d vs Group %d: p = %.4g\n', i, j, pval);
    end
end
fprintf('-------------------------------\n');
end