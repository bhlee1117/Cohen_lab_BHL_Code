classdef nidaq < dynamicprops
    properties % (Access = private)
        outChannelsConfigData
        inChannelsConfigData
        rateConfigData
        clockConfigData
        numOutChannels
        numInChannels
        syncSession
        clock
        dataIn
    end
    properties (SetAccess = private)
        lastStaticValues
    end
    methods (Access = private)
        function obj = nidaq(varargin)
            obj.outChannelsConfigData = varargin{1}; % cell string table with all channels names and configuration arguments
            obj.inChannelsConfigData = varargin{2}; % cell string table with all channels names and configuration arguments
            obj.rateConfigData = varargin{3}; % this is just a number with the approx rate of external clock input to be used
            obj.clockConfigData = varargin{4}; % cell string array of clock configuration arguments
            obj.numOutChannels = size(obj.outChannelsConfigData,1);
            obj.numInChannels = size(obj.inChannelsConfigData,1);
            obj.lastStaticValues = zeros(1,obj.numOutChannels); % initialize values to output value 0
            for it = 1:obj.numOutChannels
                propName = obj.outChannelsConfigData{it,1};
                p = obj.addprop(propName);
                p.SetObservable = true;
                addlistener(obj,propName,'PostSet',@(src,evnt)propPostSetFcn(obj,src,evnt));
            end
            obj.syncInit
            for it = 1:obj.numOutChannels
                propName = obj.outChannelsConfigData{it,1};
                obj.(propName) = 0;
            end
        end
        function syncInit(obj)
            obj.syncSession = daq.createSession('NI');
            for it = 1:obj.numOutChannels
                switch obj.outChannelsConfigData{it,2}{2}(1:2)
                    case 'Po' % for digital ports
                        obj.syncSession.addDigitalChannel(obj.outChannelsConfigData{it,2}{:});
                    case 'ao' % analog output
                        obj.syncSession.addAnalogOutputChannel(obj.outChannelsConfigData{it,2}{:});
                    otherwise
                        disp(obj.outChannelsConfigData{it,2})
                        warning 'not implemented'
                end
            end
            for it = 1:obj.numInChannels
                switch obj.inChannelsConfigData{it,2}{2}(1:2)
                    case 'Po'
                        obj.syncSession.addDigitalChannel(obj.inChannelsConfigData{it,2}{:});
                    case 'ai'
                        obj.syncSession.addAnalogInputChannel(obj.inChannelsConfigData{it,2}{:});
                    otherwise
                        disp(obj.inChannelsConfigData{it,2})
                        warning 'not implemented'
                end
            end
            obj.syncSession.Rate = obj.rateConfigData; % approx 102 kHz
%             obj.clock = obj.syncSession.addClockConnection(obj.clockConfigData{:}); % camera programmable out clock 102 kHz
        end
        function propPostSetFcn(obj,src,~)
            propIndex = strcmp(src.Name,obj.outChannelsConfigData(:,1));
            obj.lastStaticValues(propIndex) = obj.(src.Name);
%             if ~obj.syncSession.IsRunning
                obj.syncSession.outputSingleScan(obj.lastStaticValues);
%             end
            % disp(obj.lastStaticValues)
            pause(.005) % minimal wait to save shutters from destruction
        end
    end
    methods (Static)
        function obj = getInstance(varargin) % constructs singleton class
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                localObj = nidaq(varargin{:});
            end
            obj = localObj;
        end
    end
    methods
        function duration = syncQueueData(obj,outData) % outData must have the same number of columns as channels in the syncSession
            if obj.syncSession.DurationInSeconds > 0
                obj.syncSession.release
            end
            obj.syncSession.queueOutputData(outData + 0);
            duration = obj.syncSession.DurationInSeconds;
        end
        function syncStartBackground(obj) % will start and return immediately
            obj.syncSession.startBackground
        end
        function syncStartForeground(obj) % will start and return after output is completed
            obj.dataIn = obj.syncSession.startForeground;
        end
        function zeroAllOutputs(obj) % will start and return after output is completed
            obj.lastStaticValues = zeros(size(obj.lastStaticValues));
            obj.syncSession.outputSingleScan(obj.lastStaticValues);
        end
		
		function setRate(obj,rate)
			obj.syncSession.Rate = rate; 
        end
    end
end