function out = blueLUTmW(normalizedInput,max_mWcm2)
%
%   2018-2019 Vicente Parot
%   Cohen Lab - Harvard university
%
powerCalVmW = [
	0.05	0
	0.051	0.000161
	0.055	0.00181
	0.062	0.00538
	0.07	0.00959
	0.08	0.01482
	0.1	0.0254
	0.15	0.050
	0.2	0.0743
	0.3	0.118
	0.35	0.1362
	0.4	0.158
	0.5	0.195
	1	0.353
	2	0.594
	3	0.786
	4	0.952
	5	1.1
	6	1.261
	7	1.373
	8	1.489
	9	1.604
	10	1.72
];
% targetMaxV = 0.4243; % 100 mW/cm2
% targetMaxV = 0.2211; % 50 mW/cm2
% targetMaxV = 0.1332; % 25 mW/cm2
% targetMaxV = 0.0836; % 10 mW/cm2
% 
% targetMaxArbPower = interp1(powerCalVmW(:,1),powerCalVmW(:,2),targetMaxV,'spline',0);
% a = [ % intensity, targetMaxArbPower
%     100 .1680
%     50 .0844
%     25 0.0419
%     10 0.0167
%     ];
% x = [a(:,1) a(:,1)*0+1]\a(:,2);
% x = [ % Calibrated for 2x objective
%    0.001681552878179
%   -0.000021820615797
%     ];
% Calibrated on 2019-02-08
% a = [ % intensity, targetMaxArbPower
%     .108/acm2 interp1(powerCalVmW(:,1),powerCalVmW(:,2),.07,'spline',0)
%     .2641/acm2 interp1(powerCalVmW(:,1),powerCalVmW(:,2),.1,'spline',0)
%     .7433/acm2 interp1(powerCalVmW(:,1),powerCalVmW(:,2),.2,'spline',0)
%     1.1717/acm2 interp1(powerCalVmW(:,1),powerCalVmW(:,2),.3,'spline',0)
%     1.5627/acm2 interp1(powerCalVmW(:,1),powerCalVmW(:,2),.4,'spline',0)
%     1.925/acm2 interp1(powerCalVmW(:,1),powerCalVmW(:,2),.5,'spline',0)
%     2.736/acm2 interp1(powerCalVmW(:,1),powerCalVmW(:,2),.75,'spline',0)
%     3.455/acm2 interp1(powerCalVmW(:,1),powerCalVmW(:,2),1,'spline',0)
%     5.793/acm2 interp1(powerCalVmW(:,1),powerCalVmW(:,2),2,'spline',0)
%     ];
% x = [a(:,1) a(:,1)*0+1]\a(:,2);
x = [ % Calibrated on 2019-02-08 for 10x objective
   0.000324208964040
  -0.002315739472714
    ];
targetMaxArbPower = [max_mWcm2 1]*x;
out = interp1(powerCalVmW(:,2),powerCalVmW(:,1),normalizedInput*targetMaxArbPower,'spline',0);
