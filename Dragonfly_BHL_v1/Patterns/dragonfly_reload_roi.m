function [roirows, roicols] = dragonfly_reload_roi(hadamard_mode)
%
%   2018-2019 Vicente Parot
%   Cohen Lab - Harvard university
%
    switch hadamard_mode
        case 'rcamptopatch'
             roirows = 12*10+8:65*10+4; % central 2048 x 1440 with 10x  on 2019-02-08
             roicols = 13*10+6:87*10+4;
             % numel(roirows)*numel(roicols)*.9*.9*1e-8 % area in cm^2
             % 572 x 739 px^2
        case 'structural'
%             roirows = 1:768; % extended full fov 2017-03-29
%             roicols = 14*10+0:87*10+1;
             roirows = 32*10+0:43*10+3; % rat in vivo roi 2017-11-30 256
             roicols = 44*10+0:54*10+1;
        case 'activity'
%             roirows = 19*10+0:57*10+3; % extended central quad 2017-03-17
%             roicols = 32*10+0:69*10+5;
%             roirows = 18*10+0:55*10+3; % extended central quad 2017-12-18
%             roicols = 32*10+0:69*10+5;
            
                roirows = 18*10+0:55*10+3; % extended central half 2017-12-02
                roicols = 13*10+0:87*10+1;
%              roirows = 33*10+0:42*10+3; % rat in vivo roi 2017-11-30 256
%              roicols = 35*10+0:59*10+1;
%    
            % roirows = 33*10+0:57*10+3; % LTP stim roi 2017-04-07
            % roicols = 32*10+0:69*10+5;
        case 'voltage'
%                 roirows = 18*10+9:56*10+9; % extended central half 2017-12-02
%                 roicols = 13*10+0:87*10+1;
%                 roirows = (379-26:379+26)-0; % extended narrow strip 2017-12-21
%                 roirows = (379-60:379+60)-50; % extended narrow strip 2017-12-21
%                 roicols = 13*10+0:87*10+1;
%              roirows = 30*10+0:40*10+3; % rat in vivo roi 2017-11-30 256
%              roicols = 30*10+0:50*10+1;
%              roirows = 32*10+0:41*10+3; % auto cal 2018-08-03
%              roicols = 38*10+0:62*10+1;
             roirows = 32*10+0:41*10+3; % power cal 2018-11-01
             roicols = 36*10+0:64*10+1;
        case 'stim'
                roirows = 18*10+9:56*10+9; % extended central half 2017-12-02
                roicols = 13*10+0:87*10+1;
%              roirows = 33*10+0:42*10+3; % rat in vivo roi 2017-11-30 256
%              roicols = 35*10+0:59*10+1;
    end
end
