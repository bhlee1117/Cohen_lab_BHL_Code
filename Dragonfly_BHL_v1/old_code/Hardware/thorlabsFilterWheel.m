classdef thorlabsFilterWheel < handle% thorlabs FW102
    properties
        allOD = [ %default OD value
            0
            0.5
            1
            2
            3
            4
            ];
        nidaq
        DAQChannel = 'port0/line0';
        ChannelID
        currentPosition
    end

    methods
        function obj = thorlabsFilterWheel(nidaq) % removes backslash, but initialized as other motors
            obj.currentPosition = 1;
            obj.nidaq = nidaq;
            
            nDAQChannels = length(obj.nidaq.syncSession.Channels);
            for id = 1:nDAQChannels
                if strcmp(obj.DAQChannel,obj.nidaq.syncSession.Channels(id).ID)
                    obj.ChannelID = id;
                end
            end
        end
        
        function obj = moveTo(obj,OD) % 
            
            nextPosition = find(obj.allOD==str2num(OD));
            nSingleTrigger = nextPosition-obj.currentPosition;
            if nSingleTrigger>0
                % do nothing
            elseif nSingleTrigger <0
                nSingleTrigger = nSingleTrigger+length(obj.allOD);
            else
                return
            end
            
            for ii = 1:nSingleTrigger
                eval(['obj.nidaq.' obj.nidaq.outChannelsConfigData{obj.ChannelID,1} '= 1;']);
                pause(.1)
                eval(['obj.nidaq.' obj.nidaq.outChannelsConfigData{obj.ChannelID,1} '= 0;']);
                pause(1)
            end
            
            obj.currentPosition = nextPosition;
        end
    end
end