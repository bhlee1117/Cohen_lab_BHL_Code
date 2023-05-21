classdef AOTF_GH < handle
    properties
        session % serial session
        blue_freq = 75.8;
        blue_channel = 1;
        blue_max_power_level = 750;
        
    end
    methods
        function obj = AOTF_GH(comPort)
            obj.session = serial(comPort,'Terminator','CR');
%             obj.session.Timeout = .01;         
%             obj.session.InputBufferSize = 1; % get chars one by one to stop at '>'
            try
                fopen(obj.session);
                fclose(obj.session);
            catch
                error('%s could not be opened',comPort)
            end
            try
                fopen(obj.session);
                fprintf(obj.session,['ch' num2str(obj.blue_channel)]);
                fprintf(obj.session,['am' num2str(obj.blue_max_power_level)]);
                fprintf(obj.session,'mod');
                fclose(obj.session);
            catch 
                error('device not responding in %s',comPort)
            end
        end
%         function varargout = sendCommand(obj,command)
        function sendCommand(obj,command)
            fclose(obj.session); % lets just avoid fopen errors
            fopen(obj.session); % reopening ensures buffer is flushed
            fprintf(obj.session,command); % status
%             if nargout % only process output if requested
%                 response = {};
%                 str = '';
%                 out = char(fread(obj.session));
%                 while out~='>'
%                     switch out
%                         case 13
%                             response = [response; str];
%                             str = '';
%                         otherwise
%                             str = [str out];
%                     end
%                     out = char(fread(obj.session));
%                 end
%                 varargout{1} = response; % return cell array with mltiple char lines
%             end
            fclose(obj.session);
        end
%         function status = isOn(obj)
%             reply = obj.sendCommand('getldenable');
%             status = eval(reply{1});
%         end
%         function turnOn(obj)
%             obj.sendCommand('setldenable 1')
%         end
        function turnOff(obj)
            obj.sendCommand('off')
        end
%         function setPower(obj,mW)
%             assert(isnumeric(mW))
%             if mW < 100
%                 error('minimum laser power is 100 mW')
%             elseif mW < 200
%                 fprintf('laser output power set to less than 20%% nominal\n')
%                 obj.sendCommand(sprintf('setpower 0 %d',mW))
%             elseif mW <= 1000
%                 obj.sendCommand(sprintf('setpower 0 %d',mW))
%             else
%                 error('maximum laser power is 1000 mW')
%             end
%         end
%         function mW = getPowerSetpoint(obj)
%             reply = obj.sendCommand('getpower 0');
%             mW = eval(reply{1});
%         end
%         function mW = getPowerOutput(obj)
%             reply = obj.sendCommand('power 0');
%             mW = eval(reply{1});
%         end
    end
end