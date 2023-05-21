function ctrl_input = power_cal(mW_cm2,color,input)
% 1 mm 2x2 binning
x1 = [549 390];
x2 = [552 613];


pixel2mm = 1/sqrt(sum((x2-x1).^2));

% 480 area
width480 = 1424-810;
length480 = 1267-608;
area480 = width480*length480/4;
% 560 area

area480=0.2^2*pi;
area560 = area480;


% 480 V vs uW

ctrl_voltage_480 = xlsread('power_cal.xlsx','Sheet1','D3:D27');
power_480 = xlsread('power_cal.xlsx','Sheet1','E3:E27');
% mW_cm2_480 = power_480*1e-3/(area480*pixel2mm^2*1e-2);
mW_cm2_480 = power_480*1e-3/(area480);

% 560 mW (input) vs uW (readout)
ctrl_power_560 = xlsread('power_cal.xlsx','Sheet1','G3:G7');
power_560 = xlsread('power_cal.xlsx','Sheet1','H3:H7');
% mW_cm2_560 = power_560*1e-3/(area560*pixel2mm^2*1e-2);
mW_cm2_560 = power_560*1e-3/(area560);
%%
% figure
% plot(ctrl_voltage_480,mW_cm2_480,'b')
% % 
% figure
% plot(ctrl_power_560,mW_cm2_560,'g')
switch input
    case 'control'
        switch color
            case 'blue'
                ctrl_input = interp1(ctrl_voltage_480,mW_cm2_480,mW_cm2,'spline',0);
            case 'lime'
                ctrl_input = interp1(ctrl_power_560,mW_cm2_560,mW_cm2,'linear','extrap');
        end
    case 'intensity'
        switch color
            case 'blue'
                ctrl_input = interp1(mW_cm2_480,ctrl_voltage_480,mW_cm2,'spline',0);
            case 'lime'
                ctrl_input = interp1(mW_cm2_560,ctrl_power_560,mW_cm2,'linear','extrap');
        end
end

