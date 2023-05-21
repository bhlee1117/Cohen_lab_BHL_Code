classdef APTmotor < handle
    properties
        h
    end
    methods
        function obj = APTmotor(serial,position,hFig)
            if ~exist('hFig','var')
                hFig = prepareFigure;
            end
            if ~exist('position','var')
                position = [0 0,500,300];
            end
            obj.h = actxcontrol('MGMOTOR.MGMotorCtrl.1',position,hFig);
            obj.h.HWSerialNum = serial;
            obj.h.StartCtrl;
        end
        function currentPosition = position(obj)
            currentPosition = obj.h.GetPosition_Position(0);
        end
        function moveTo(obj,nextPosition,timeout)
            if ~exist('timeout','var')
                timeout = 0;
            end
            assert(isnumeric(nextPosition),'moveTo position must be numeric')
            [~, oldMinLim, oldMaxLim, ~, ~, ~] = obj.h.GetStageAxisInfo(0,0,0,0,0,0);
            if nextPosition < oldMinLim
                error('target position below lower software limit') % better stop here than in APT
            end
            if nextPosition > oldMaxLim
                error('target position above higher software limit') % better stop here than in APT
            end
            distance = (nextPosition - obj.position); % * obj.distancePositionRatio;
            obj.h.SetRelMoveDist(0,distance);
            obj.h.MoveRelative(0,0);
            startTic = tic;
            elapsedTime = toc(startTic);
            while elapsedTime < timeout
                if ~obj.isMoving
                    break
                end
                elapsedTime = toc(startTic);
            end            
        end
        function movingFlag = isMoving(obj)
            val = obj.h.GetStatusBits_Bits(0);
            movingFlag = ~~bitget(typecast(val,'uint64'),33);
            % % debug
%             fprintf('%d',bitget(typecast(val,'uint64'), 1:16)), fprintf(' ')
%             fprintf('%d',bitget(typecast(val,'uint64'),17:32)), fprintf(' ')
%             fprintf('%d',bitget(typecast(val,'uint64'),33:48)), fprintf(' ')
%             fprintf('%d',bitget(typecast(val,'uint64'),49:64)), fprintf('\n')
        end
    end
    methods (Static)
        function hFig = prepareFigure(position)
            hFig = figure('Visible','off');
            hFig.NumberTitle = 'off';
            hFig.Name = 'APT motor';
            hFig.ToolBar = 'none';
            hFig.MenuBar = 'none';
            if exist('position','var')
                hFig.InnerPosition = position;
            else
                hFig.InnerPosition = hFig.InnerPosition.*[1 1 0 0] + [0 0,500,300];
            end
            hFig.Resize = 'off';
            hFig.HandleVisibility = 'off';
            hFig.Visible = 'on';
        end
    end
end