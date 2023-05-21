classdef slm_device


    properties
            %Device Characteristics
    end
    methods(Static)
        function obj=slm_device() %
            obj.init
        end
        
        function init(obj)
            dll_path = 'C:\Program Files\Meadowlark Optics\Blink 1920 HDMI\SDK\Blink_C_wrapper.dll';
            h_path = 'C:\Program Files\Meadowlark Optics\Blink 1920 HDMI\SDK\Blink_C_wrapper.h';
            if libisloaded('Blink_C_wrapper')
                fprintf('Library was loaded...\nunload now...')
                unloadlibrary('Blink_C_wrapper');
            end
            loadlibrary(dll_path,h_path);

            % Matlab automatically escapes backslashes (unlike most languages)
            % Linear is a default LUT, you should replace this path to your customer
            % LUT that was delivered with your SLM. This calibrates the nonlinear
            % response of the LC to voltage.
            lut_file = 'C:\Program Files\Meadowlark Optics\Blink 1920 HDMI\LUT Files\linear.lut';

            calllib('Blink_C_wrapper', 'Create_SDK');

            disp('Blink SDK was successfully constructed');

            % Load the lookup table to the controller. Make sure that your com port is
            % correct. You can open Blink to validate the COM port. 
            calllib('Blink_C_wrapper', 'Load_lut', lut_file);
            calllib('Blink_C_wrapper','Set_channel',0);
            
        end
        
        function free(obj)
           calllib('Blink_C_wrapper', 'Delete_SDK');
           unloadlibrary('Blink_C_wrapper');
        end
    end
    
    methods
        function project(obj,pat)
            pause(.1) % wait time:free up GPU for SLM
            if ~isa(pat,'uint8')
                pat = uint8(mat2gray(pat)*255);
            end
            
            calllib('Blink_C_wrapper', 'Write_image', rot90(pat,3), true);
            disp('project success!')

        end
    end
    
end