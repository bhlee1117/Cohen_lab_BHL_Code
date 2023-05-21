function chr = num2ordinal(num,sfx)
% Convert a numeric array to a character array of integers and ordinal suffixes.
%
% (c) 2014-2020 Stephen Cobeldick
%
% Convert a numeric array into a character array containing integers and
% ordinal suffixes, e.g. 1 -> '1st', rounding any fractional values to the
% nearest integers. NUM2ORDINAL returns the correct suffixes for values
% ending in 11, 12, or 13 (unlike the Image Toolbox's IPTNUM2ORDINAL).
%
%%% Syntax:
% chr = num2ordinal(nun)
% chr = num2ordinal(num,sfx)
%
%% Examples %%
%
% >> num2ordinal(1)
% ans = '1st'
%
% >> num2ordinal(1:6)
% ans = ['1st';'2nd';'3rd';'4th';'5th';'6th']
%
% >> num2ordinal([1,11,111,1111])
% ans = ['   1st';'  11th';' 111th';'1111th']
%
% >> num2ordinal(100:113,true)
% ans = ['th';'st';'nd';'rd';'th';'th';'th';'th';'th';'th';'th';'th';'th';'th']
%
% >> num2ordinal(intmax('int64')-4)
% ans = '9223372036854775803rd'
%
% >> vec = [-1,-0,0;-Inf,NaN,Inf];
% >> reshape(strtrim(cellstr(num2ordinal(vec))),size(vec))
% ans = {'-1th','-0th','0th'
%      '-Infth','NaNth','Infth'}
%
%% Input & Output Arguments %%
%
%%% Inputs (*==default):
% num = NumericArray (any size), with values to convert to string.
% sfx = LogicalScalar, true/false* -> suffix only / integer & suffix.
%
%%% Output:
% chr = Character Array, with integers of rounded <num> values with
%       ordinal suffixes. The rows are linear indexed from <num>.
%
% See also NUM2WORDS WORDS2NUM NATSORT NUM2BIP NUM2SIP SPRINTF INT2STR NUM2STR NUM2CELL CELLSTR STRTRIM STRCAT INTMAX

%% Input Wrangling %%
%
typ = class(num);
switch typ(1:3)
    case {'dou','sin'}
        assert(isreal(num),'SC:num2ordinal:NotRealNumeric',...
			'First input <num> must be a real numeric.')
        fmt = '.0f';
    case 'uin'
        fmt = 'lu';
    case 'int'
        fmt = 'ld';
    otherwise
        error('SC:num2ordinal:NotNumericInput',...
			'First input <num> must be a numeric scalar/vector/matrix/array.')
end
num = num(:);
%
%% Convert from number to character array plus 'th' suffix %%
%
if isempty(num) % Empty Numeric:
    chr = '';
    return
    %
elseif nargin>1 && sfx % Suffixes Only:
    two = 2;
    chr(1:numel(num),2) = 'h';
    chr(:,1) = 't';
    %
elseif isscalar(num) % Numeric scalar:
    chr = sprintf(['%',fmt,'th'], num);
    two = numel(chr);
    %
else % Numeric Matrix:
    % Calculate the number of digits required:
    two = 1 + max(0, floor(log10(double(abs(num)))));
    two(~isfinite(num)) = 3;
    two = 2 + max(two + (num<0 | 1./num<0));
    % Convert from numeric to character:
    fmt = sprintf('%%%.0f%sth', two-2, fmt);
    chr = reshape(sprintf(fmt, num), two, []).';
    %
end
%
%% Determine Ordinal Suffixes %%
%
rm = rem(num,100);
th = 10.5<=rm & rm<13.5;
rm = rem(num,10);
st = ~th & 0.5<=rm & rm<1.5;
nd = ~th & 1.5<=rm & rm<2.5;
rd = ~th & 2.5<=rm & rm<3.5;
%
% Replace 'th' suffixes with the correct suffixes:
one = two-1;
chr(st,one) = 's';
chr(st,two) = 't';
chr(nd,one) = 'n';
chr(nd,two) = 'd';
chr(rd,one) = 'r';
chr(rd,two) = 'd';
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%num2ordinal