classdef flexSineWave < baseWave
    
    %flexSineWave Flexible sine waveform. Adds detailed amplitude control
    %   to base waveforms.
    %
    %   Waveforms are represented by any combination of defining
    %   parameters. All derived properties can be accessed without storing
    %   underlying waveform matrices. This allows for efficient handling of
    %   large waveforms. Waveforms can be processed with standard math
    %   operations and concatenated. Calibrated waveforms can be defined
    %   using a look-up table.
    %
    %   Methods for class flexSineWave:
    %     flexSineWave
    %       plus all methods inherited from baseWave
    %
    %   Properties for class flexRampWave:
    %     offset
    %     initialHeight
    %     finalHeight
    %     baseAmp
    %     initialAmp
    %     finalAmp
    %       plus all properties inherited from baseWave
    %
    %   See also baseWave, rampWave, flexRampTrainWave, arbWave.
    %
    %   2018-2021 Vicente Parot, Yitong Qi
    %   Cohen Lab - Harvard university
    
    properties
        phase = 0
        tBefore = 0;
        tOn = 0;
        offset = 0;
        period = 0;
        wAmplitude = 1;
    end
    
    methods
        function obj = flexSineWave(varargin)
            switch nargin
                case 0
                case 1
                    obj.dt = varargin{1};
                case 2
                    obj.dt = varargin{1};
                    obj.duration = varargin{2};
            end
        end
        
        
        function value = localGetAmplitude(obj)
            if ...
                    ~isempty(obj.tOn) && ...
                    ~isempty(obj.tBefore) && ...
                    ~isempty(obj.phase) && ...
                    ~isempty(obj.offset)
                value = sin(obj.time/obj.period*2*pi)*obj.wAmplitude;
                value = value .* (obj.tBefore < obj.time & obj.time <= (obj.tBefore+obj.tOn));
                value = value + obj.offset;
            else
                error('waveform paremters underdetermined')
            end
        end
    end
    
end