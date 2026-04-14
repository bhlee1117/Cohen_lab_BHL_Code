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

p = inputParser;
p.addParameter('Method','auto',@(s)ischar(s)||isstring(s));
p.addParameter('TestFcn',[],@(f)isempty(f)||isa(f,'function_handle'));
p.parse(varargin{:});
method  = lower(string(p.Results.Method));
testFcn = p.Results.TestFcn;

if isPaired
    numGroups = size(data, 2);
else
    numGroups = numel(data);
end

p_values = nan(numGroups);

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

        % ----- run requested test -----
        if ~isempty(testFcn)
            pval = testFcn(x,y);

        else
            if method == "auto"
                method = isPaired * "ttest" + ~isPaired * "ranksum"; % default behavior
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
    end
end
end