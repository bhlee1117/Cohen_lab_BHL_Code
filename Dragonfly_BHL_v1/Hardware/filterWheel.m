classdef filterWheel < APTmotor
    properties
        allFilters = { % one ring to rule them all
            'Arch'
            'GFP'
            'RCaMP'
            'RFP'
            'CFP'
            };
        motorAmplitudes = [ % calibrated by adjusting after stepping multiple turns
            -2/(567/570.54)/5*.8*(570.54-4.4)/570.54
            -1/(567/570.54)/5*.8*(570.54-3.6)/570.54
             0
            +1/(567/570.54)/5*.8*(570.54-7.5)/570.54
            +2/(567/570.54)/5*.8*(570.54-6.0)/570.54
            ];
        currentFilter
    end
    methods
        function obj = filterWheel(varargin) % removes backslash, but initialized as other motors
            obj = obj@APTmotor(varargin{:});
            obj.h.SetBLashDist(0,0);
            obj.h.GetBLashDist_BLashDist(0);
        end
        function moveTo(obj,nextPosition,varargin) % accepts a position in mm, optionally a timeout value 
            % set filter wheel limits above and below current target
            assert(isnumeric(nextPosition),'moveTo position must be numeric. try switchTo?')
            possibleValues = [nextPosition obj.position 0];
            newMinLim = min(possibleValues)-1;
            newMaxLim = max(possibleValues)+1;
            [a, ~, ~, units, pitch, f] = obj.h.GetStageAxisInfo(0,0,0,0,0,0);
            obj.h.SetStageAxisInfo(a, newMinLim, newMaxLim, units, pitch, f);
            obj.moveTo@APTmotor(nextPosition,varargin{:})
        end
        function switchTo(obj,nextFilter,varargin) % accepts a position string, optionally a timeout value
            if isempty(obj.currentFilter)
                obj.currentFilter = obj.allFilters{3};
            end
            originIndex = strcmp(obj.allFilters,obj.currentFilter);
            destinationIndex = strcmp(obj.allFilters,nextFilter);
            assert(sum(originIndex)==1) % there must be only one index to find them
            assert(sum(destinationIndex)==1)
            numberOfPositions = numel(obj.allFilters);
            forwardCount = 1:+1:numberOfPositions;
            bwkwardCount = numberOfPositions:-1:1;
            forwardSteps = mod(forwardCount-forwardCount(destinationIndex),numberOfPositions);
            bwkwardSteps = mod(bwkwardCount-bwkwardCount(destinationIndex),numberOfPositions);
            stepsToDestination = forwardSteps.*(forwardSteps<bwkwardSteps)-bwkwardSteps.*(forwardSteps>bwkwardSteps);
            motorAmplitudeIndex = stepsToDestination(originIndex) + 3;
            obj.moveTo(obj.position + obj.motorAmplitudes(motorAmplitudeIndex),varargin{:})
            obj.currentFilter = obj.allFilters{destinationIndex};
        end
        function resetCurrentFilterPositionTo(obj,trueCurrentFilter) % resets the internal state of current filter
            if ~exist('trueCurrentFilter','var')
                trueCurrentFilter = 'GFP';
            end
            filterIndex = strcmp(obj.allFilters,trueCurrentFilter);
            assert(sum(filterIndex)==1) % there must be only one filter to bring them all and in the darkness bind them
            obj.currentFilter = obj.allFilters{filterIndex};
        end
    end
end