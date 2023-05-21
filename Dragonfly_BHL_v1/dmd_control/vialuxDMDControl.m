classdef vialuxDMDControl < handle
    properties
        api
        device
        seq
        roirows
        roicols
        type
    end
    methods (Access = private)
        function obj = vialuxDMDControl(varargin)
            obj.init
        end
    end
    methods (Static)
        function obj = getInstance(varargin)
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                localObj = vialuxDMDControl(varargin{:});
            end
            obj = localObj;
        end
    end
    methods
        function init(obj)
            %% initialize DMD

            % 2016 Vicente Parot
            % 2021 Yitong Qi
            % Cohen Lab - Harvard University

            % Check dmd api version
            if ~exist('api_version.txt','file')
                fig = uifigure;
                alp_version = uiconfirm(fig,'Please check api version','','option',...
                {'alpV42x64',...
                'alpV43x64'},...
                'closefcn',@(h,e) close(fig));
            
                fid = fopen('api_version.txt','w');
                fwrite(fid,alp_version)
                fclose(fid);
            else
                alp_version = fileread('api_version.txt');
            end
            switch alp_version
                case 'alpV42x64'
                    if isempty(obj.api) || ~isa(obj.api,'alpV42x64')
                        obj.api = alpload('alpV42x64');
                    end
                case 'alpV43x64'
                    if isempty(obj.api) || ~isa(obj.api,'alpV43x64')
                        obj.api = alpload('alpV43x64');
                    end
                otherwise
                    error 'unknown rig'
            end
            % connects to device, resets connection if there already was one
            obj.device = alpdevice(obj.api);
            switch obj.device.alloc
                case obj.api.DEFAULT
                    disp 'device alloc ok'
                case obj.api.NOT_ONLINE
                    disp 'device not online'
                otherwise
                    display(['device alloc returned ' num2str(alloc_val)])
            end
        end
        function setType(obj,type)
            obj.type = type;
        end
        function free(obj)
            obj.device.free
        end
        function reload_roi(obj,type)
%             clear functions
            [obj.roirows, obj.roicols] = dragonfly_reload_roi(type);
        end
%         function load_patterns(obj)
%             obj.reload_roi(obj.type)
%             
% %             clear functions
%             % new_labview_patterns_fingerprint = [hash_fcn(@dragonfly_generate_alp_patterns) hash_fcn(@dragonfly_reload_roi)];
%             % if ~exist('labview_patterns_fingerprint') || ~isequal(labview_patterns_fingerprint,new_labview_patterns_fingerprint)
%                     alp_patterns = dragonflyGenerateAlpPatterns(obj.device,obj.roirows,obj.roicols,obj.type);
%             %        labview_patterns_fingerprint = [hash_fcn(@dragonfly_generate_alp_patterns) hash_fcn(@dragonfly_reload_roi)];
%             % end
% %             dragonfly_load_sequence        
%             obj.load_sequence(alp_patterns)
%         end
        function load_sequence(obj,alp_patterns)
%             clc
            obj.device.free;
%             disp 'device free' 
            obj.init
%             lab_init_device
            % ExposureMilliseconds = 10;
            % CameraLines = 2048;
            ModeStringSyncOrUninterrupted = 'uninterrupted';
            % setup_load_start_sequence
            %%
            obj.device.stop;
            obj.device.halt;

            BitPlanes = 1;
            PicOffset = 0;
            PicNum = size(alp_patterns,3);

            % PictureTimeSpacing = ceil(ExposureMilliseconds*128/5); % [us]
            TriggerInDelay = obj.api.DEFAULT; % -30 + ceil(9.75*ceil(CameraLines/2)); % [us] only relevant in slave mode
%             PictureTimeExcess = 25; % [us] % must be > 2 us per instructions. not a dmd parameter, for parameter calculation only

            TriggerSynchDelay = obj.api.DEFAULT; % [us] only relevant in master mode
            TriggerSynchPulseWidth = obj.api.DEFAULT; % [us]
            PictureTime = 44; % api.DEFAULT; % floor(ExposureMilliseconds*1000/1.023)-PictureTimeSpacing; % [us]
            IlluminateTime = obj.api.DEFAULT; % PictureTime-TriggerInDelay-PictureTimeExcess-TriggerSynchDelay; % [us]
%             fprintf('PictureTime: %d us\nIlluminateTime: %d us\n',PictureTime,IlluminateTime)

            obj.seq = alpsequence(obj.device);
            obj.seq.alloc(BitPlanes,PicNum);
            obj.seq.control(obj.api.DATA_FORMAT,obj.api.DATA_BINARY_TOPDOWN);
            switch ModeStringSyncOrUninterrupted
                case 'uninterrupted'
                    obj.seq.control(obj.api.BIN_MODE,obj.api.BIN_UNINTERRUPTED); % to display the pattern
                    % until next trigger, even in slave mode, regardless of IlluminateTime.
                    % Watch for bleedthrough into the rolling shutter due to delayed
                    % responsivity.  
                case 'sync' % do nothing
                otherwise % do nothing
            end
            obj.seq.timing(IlluminateTime, PictureTime, TriggerSynchDelay, TriggerSynchPulseWidth, TriggerInDelay);
            [~, PictureTime] = obj.seq.inquire(obj.api.PICTURE_TIME);
            [~, IlluminateTime] = obj.seq.inquire(obj.api.ILLUMINATE_TIME);
            fprintf('PictureTime: %d us\nIlluminateTime: %d us\n',PictureTime,IlluminateTime)
            fprintf('loading %d patterns ... \n',PicNum)
            obj.seq.put(PicOffset,PicNum,permute(alp_patterns,[3 1 2]));

            obj.device.projcontrol(obj.api.PROJ_MODE,obj.api.SLAVE_VD);
            obj.device.startcont(obj.seq);
            fprintf([8 'done\n'])

        end
        function all_on(obj)
            pat = ones(obj.device.height,obj.device.width);
            obj.project(pat);
        end
        function all_off(obj)
            pat = zeros(obj.device.height,obj.device.width);
            obj.project(pat);
        end
        function dots(obj)
            pat = rand(obj.device.height,obj.device.width)<.05;
            obj.project(pat);
        end
        function special(obj)
            pat = rand(obj.device.height,obj.device.width)<.5;
            obj.project(pat);
        end
        function roi(obj)
            pat = zeros(obj.device.height,obj.device.width);
            pat(obj.roicols,obj.roirows) = 1;
            obj.project(pat);
        end
        function project(obj,binaryPattern)
            try
            obj.device.stop;
            obj.device.halt;
            obj.seq.free;
            end
            PicNum = int32(1);
            PicOffset = int32(0);
            BitPlanes = int32(1);
            obj.seq = alpsequence(obj.device);
            obj.seq.alloc(BitPlanes,PicNum);
            % seq.control(api.DATA_FORMAT,api.DATA_BINARY_TOPDOWN);
            obj.seq.control(obj.api.BIN_MODE,obj.api.BIN_UNINTERRUPTED);
            obj.seq.timing(10E6-2E-6, 10E6, 0, 0, 0);
            obj.seq.put(PicOffset, PicNum, binaryPattern*255);
            obj.device.projcontrol(obj.api.PROJ_MODE,obj.api.MASTER);
            obj.device.startcont(obj.seq);
%             obj.device.put(binaryPattern*255)
%             [alp_return, projection_mode] = obj.device.projinquire(obj.api.PROJ_MODE);
%             if(projection_mode ~= obj.api.MASTER)
%                     obj.device.halt;
%                     obj.device.projcontrol(obj.api.PROJ_MODE,obj.api.MASTER);
%                     obj.device.put(binaryPattern*255)
%             end                    
        end
        function result = status(obj)
            result = ~isempty(obj.api) && ~isempty(obj.device) && ~obj.device.alloc;
        end
    end
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           