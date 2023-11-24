function publ(width)
%
% change default linewidths of all new and present figures
% to make publication quality figures
%
% call publ(width)
%      default value for width is 0.5
%      typical value for width for figure to be used in a publication is 1-2

%
% Sverre Holm, University of Oslo 6. Nov 1998
% Modified by Adam Cohen 9 March 2009 to work with complex graphics

if (nargin ==1),
% change all following figures: 	
	set(0,'DefaultAxesLineWidth',width);
	set(0,'DefaultLineLineWidth',width);
	set(0,'DefaultPatchLineWidth',width);
	get(0,'Default');
%
% change this figure:
	set(gca,'LineWidth', width);
    
    set(gcf,'DefaultLineLineWidth',width);
    %	h = get(gca,'children');
%   set(h,'LineWidth',width)
end
return
