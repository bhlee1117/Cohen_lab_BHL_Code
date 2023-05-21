function mWcm2 = blueLUT_V_2_mWcm2(v)
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
arbPower = interp1(powerCalVmW(:,1),powerCalVmW(:,2),v,'spline',0);
x = [ % Calibrated on 2019-02-08 for 10x objective
   0.000324208964040
  -0.002315739472714
    ];
mWcm2 = arbPower/x(1);