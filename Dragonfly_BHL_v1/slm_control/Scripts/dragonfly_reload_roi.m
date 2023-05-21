function [roirows, roicols] = dragonfly_reload_roi(hadamard_mode)
    switch hadamard_mode
        case 'structural'
            roirows = 1:768; % extended full fov 2017-03-29
            roicols = 14*10+0:87*10+1;
%              roirows = 33*10+0:42*10+3; % rat in vivo roi 2017-11-30 256
%              roicols = 35*10+0:59*10+1;
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
             roirows = 33*10+0:42*10+3; % rat in vivo roi 2017-11-30 256
             roicols = 35*10+0:59*10+1;
        case 'stim'
                roirows = 18*10+9:56*10+9; % extended central half 2017-12-02
                roicols = 13*10+0:87*10+1;
%              roirows = 33*10+0:42*10+3; % rat in vivo roi 2017-11-30 256
%              roicols = 35*10+0:59*10+1;
    end
end
