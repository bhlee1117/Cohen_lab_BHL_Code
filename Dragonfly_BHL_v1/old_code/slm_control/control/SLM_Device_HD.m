classdef SLM_Device %May need to add option arguments to classdef later
%Cohen Lab SLM_Device Controller (prototype version)
%Authors: Hunter Davis and David Wong-Campos
%Contact: huntercoledavis@gmail.com,jwongcampos@gmail.com

%This is a polymorhpic SLM_Device Controller containing all necessary
%methods and data structures for generating holograms from a variety of
%target types including point-based and region-based matlab ROIs.
%Additional methods are available for calibrating the transformation matrix
%between a camera object and the SLM. As the confocal scanners will
%be registered as devices with their own transforms to the camera, I've
%added an additional feature that allows the generation of SLM holograms
%based on images prodcued from a confocal scan.

%Remaining critical tasks include but are not limited to:
%Including wavefront correction in the SLM hologram generator, replacing
%the enforced number of iterations on the hologram generator with an error
%tolerance and a max number of iterations,attaching an NI DAQ device
%digital output as a hardware trigger(not available on all SLMs, so there
%should be an object property specifying if this is allowed), writing some
%example scripts that make use of this controller, and writing the
%mid-level wrapper that will integrate this controller with the controllers
%for other devices on a rig.

%Cool features that we should add in the future include the option to generate 3D holograms
%based on a stack of ROIs, which could be produced from an input confocal
%stack and the option to generate exotic light with optical vortices etc.

%None of the filepaths have been fully fleshed out as we still need to decide
%how we are going to structure the final file tree of the repository.

    properties
            %Device Characteristics
            Device_Model
            Dimensions
            LUT_Filepath

            %Hologram Information
            Blaze
            WavefrontCorrection
            FocusX
            FocusY
            Z_GS
            Hologram_Stack
            Stack_Frame_Number
            Spot_Width
            Target

            %Transform Calibration Information
            tform_file
            calpoints_file
            tform
            calpoints
            calPS
            calthresh
            
            demo_mode
     end

    methods
        function obj=SLM_Device(basepath,demo_mode) %
              %Start by loading the parameters file, which contains the
              %basic information about the SLM device in addition to
              %filepaths that direct the driver to data required for
              %initialization.
              params=string(readcell(fullfile(basepath,'Rig_Params','SLM.txt'),'Delimiter',':'));
              obj.Device_Model=params(params(:,1)=='SLM_Device',2);
              obj.Dimensions=str2num(params(params(:,1)=='Dimensions',2));
              obj.Blaze=str2num(params(params(:,1)=='Default_Blaze',2));
              obj.LUT_Filepath=params(params(:,1)=='LUT_Filepath',2);
              obj.Spot_Width=str2num(params(params(:,1)=='Spot_Width',2)); %approximate gaussian width of spot on SLM face
              obj.tform_file=params(params(:,1)=='tform_file',2);
              obj.tform=affine2d(cell2mat(readcell(obj.tform_file)));
              obj.calpoints_file=params(params(:,1)=='calpoints_file',2);
              obj.calPS=params(params(:,1)=='calPS',2);
              obj.calthresh=10;
              obj.demo_mode=demo_mode;
              %Now load the calibration points and expand them to the SLM
              %dimensions.
              fractional_calpoints=cell2mat(readcell(obj.calpoints_file));
              obj.calpoints=ceil(fractional_calpoints.*obj.Dimensions);

              %A new functionality with this driver is the conveinient
              %storage of multiple pre-calculated holograms inside of the
              %SLM Device object using a structure called "Hologram_Stack".
              %See the methods list for how to generate and append new
              %holograms as well as iterate through the stack. The
              %"Hologram_Stack" is externally callable and visible, so a
              %parent program can directly write to it or read from it as
              %needed.

              obj.Hologram_Stack=struct();
              obj.Stack_Frame_Number=1;
              obj.Target=zeros(obj.Dimensions);
              %Unfortunately, the matlab driver for Meadowlark SLMs is
              %not very flexible due to Matlab's inability to handle
              %pointer structures. We don't get a useable handle out of the
              %SDK intialization. This isn't a problem for single SLM rigs,
              %but if we ever try to work on a multi-SLM
              %rig, we may need to customize the Meadowlark Driver.
              %Initialize SDK loads all the necessary dlls, libraries, and
              %LUTs.
              if obj.demo_mode==0
                InitalizeSDK(obj.LUT_Filepath);
              end
              obj.WavefrontCorrection=[]; %Need to add WFC to SLM_routine
              obj.FocusX=[];
              obj.FocusY=[];
              

        end
        %SLM_genholo is the base method for generating holograms from a
        %given target. The only required inputs are the SLM_Device object
        %(obj) and a Target image. Additional options can be modified using
        %a namae value pair. (e.g. SLM_genholo(obj,Target,'Blaze',[1 1])
        %would impose a blaze of [1 1] for the hologram.) Default values
        %for all options are specified in the arguments block. Instead of
        %automatically writing to the SLM. The user has the option to store
        %the hologram as a frame in the Hologram_Stack by setting the
        %"append_to_stack" option to 1. By default, the stack is cleared
        %when a new hologram is written. This ensures that for the "simple"
        %standard use case of projecting a single image for an experiment,
        %the user is sheilded from needing to worry about iterating through
        %the stack. Setting "write_when_complete" to 1 or using the
        %"Display_Simple" method are probably the best ways to exploit the
        %controller for the simple use case.
        function obj=SLM_genholo(obj,Target,options)
          arguments%The "arguments" block is new to matlab for 2019b. Please see mathworks.com for more info.
            obj SLM_Device
            Target (:,:) double
            options.Blaze (1,2) {mustBeNumeric} = obj.Blaze
            options.Focus (1,1) {mustBeNumeric} = 0
            options.gpu_avail (1,1) logical = numel(gpuDevice)>0
            options.Spot_Width (1,1) {mustBeNumeric} = obj.Spot_Width
            options.numiterations   (1,1) {mustBeNumeric} = 50
            options.write_when_complete (1,1) logical = 0
            options.append_to_stack (1,1) logical = 0
            options.clear_stack (1,1) logical = 1
            options.progress_bar_enable logical=0
            options.progress_bar matlab.ui.control.LinearGauge
          end
          obj.Target=Target;
          if options.progress_bar_enable==1
            obj.Z_GS=SLM_routine(Target,options.Spot_Width,options.numiterations,options.gpu_avail,'progress_bar',options.progress_bar);
          else
            obj.Z_GS=SLM_routine(Target,options.Spot_Width,options.numiterations,options.gpu_avail);
          end
          disp(size(obj.Z_GS))
          if options.write_when_complete
            obj.SLM_Project();
          end
          
          if options.append_to_stack
            if isempty(obj.Hologram_Stack)
                obj.Hologram_Stack(1).frame=obj.Z_GS;
            else
                obj.Hologram_Stack(numel(obj.Hologram_Stack)+1).frame=obj.Z_GS;
            end
          end

          if options.clear_stack
            obj.Hologram_Stack=struct();
            obj.Hologram_Stack(1).frame=obj.Z_GS;
            obj.Stack_Frame_Number=1;
          end
        end

        %DisplayNextFrame will project the hologram at the Stack_Frame_Num
        %position in the hologram stack and advance the Stack_Frame_Num
        %value by 1.

        function obj=DisplayNextFrame(obj)
            SLM_Project(obj.Hologram_Stack(obj.Stack_Frame_Number));
            obj.Stack_Frame_Number=mod(obj.Stack_Frame_Number,...
                numel(obj.Hologram_Stack))+1;
        end
        %Display_Simple projects the first hologram in the stack to the SLM.
        %For many uses of the controller, the stack will only contain one
        %hollogram at a time and Display_Simple will be the main method of
        %projecting the hologram on the SLM.
        function obj=Display_Simple(obj)
            SLM_Project(obj.Hologram_Stack(1));
        end

        %SLM_TCal calibrates the affine transformation that registers the camera
        %coordinates to equivalent SLM coordinates. Requires a "Camera_Device"
        %(Cam) as an input, which is a class we have yet to write. First, the
        %controller projects the control points using the SLM. Then, the camera
        %will snap an image and pass it to the SLM device in the array %"testimage". Then, the Target in SLM coordinates ("fixed") is registered
        %with the snapped camera image using the Transform_Cal function. See
        %regirstation_guide.md in docs for explanation of this function. The
        %returned 2daffine transform object "tform_out" is passed to the SLM
        %device object and saved to the tform file.

        function obj=SLM_TCal(obj,Cam)
          arguments
           obj  SLM_Device
           Cam  Camera_Device
          end
            fixed=obj.GenTargetFromSpots(obj.calpoints,obj.Dimensions,1);
            [obj,holo]=SLM_genholo(obj,fixed);
            SLM_Project(holo);
            testimage=Cam.Snap();
            tform_out=Transform_Cal(fixed,testimage,obj.calPS,obj.calthresh);
            obj.tform=tform_out;
            fid=fopen(obj.tform_file,'wt');
            for ii = 1:size(tform_out.T,1)
                fprintf(fid,'%16f\t',tform_out.T(ii,:));
                fprintf(fid,'\n');
            end
            fclose(fid);
        end


        function obj=SLM_ManCal(obj,tax,testimage)
          arguments
           obj         SLM_Device
           tax         
           testimage
          end
            fixed=obj.GenTargetFromSpots(obj.calpoints,obj.Dimensions,1);
            figure
            imagesc(tax,testimage);
            roi=drawpolyline(tax);
            xpts=roi.Position(:,1);
            ypts=roi.Position(:,2);
            t_est = estimateGeometricTransform([xpts ypts],[obj.calpoints(:,2) obj.calpoints(:,1)],'affine');
            tform_out=Transform_Cal(fixed,testimage.testim,20,3,1,t_est);
            obj.tform=tform_out;
            fid=fopen(obj.tform_file,'wt');
            for ii = 1:size(tform_out.T,1)
                fprintf(fid,'%16f\t',tform_out.T(ii,:));
                fprintf(fid,'\n');
            end
            fclose(fid);
        end
%GenROI allows the user to create ROI targets in SLM coordinates from images
%produced by either a camera or a confocal scan. See registration_guide.md for an
%overview of how the SLM to confocal scan registration works. The  "smoothing"
%option sets the width of the gaussian blur applied to the Target in SLM units. A %finite blur greatly improves the stability of the hologram generator without
%seriously degrading the resolution. The "ROI_Type" option can be set to
%'points', 'polyedge', or 'polyregion'. If ROI_Type type is set to 'points', then
%the user will be prompted to click on each position on the image where they
%would like to drop a spot. If ROI_Type is set to 'polyedge' the user creates an
%unfilled polygonal ROI. This could be useful for selectively stimulating the
%cell membrane. Finally, the 'polyregion' ROI_Type allows the user to generate
%filled ROIs. The "append" option allows the user to add the generated target
%onto the stored object target, which means you can create a target from mixed
%ROI types. The renormalization of the targets after the addition will be
%relevant if there is overlap between the generated target and the stored target.

        function [obj,Target]=GenROI(obj,ax,testimage,options)
          arguments
            obj SLM_Device
            ax
            testimage
            options.smoothing (1,1) double = 50
            options.ROI_Type string = 'points'
            options.append logical = 0
          end
                transform=obj.tform;
                I=imagesc(ax,testimage);
                ax.XLim = [0 I.XData(2)];
                ax.YLim = [0 I.YData(2)];
           switch options.ROI_Type
             case 'points'
                roi=drawpolyline(ax);
                i=roi.Position(:,1);
                j=roi.Position(:,2);
                [x,y]=transformPointsForward(transform,i,j);
                Target=obj.GenTargetFromSpots([y,x],obj.Dimensions,options.smoothing);
             case 'polyedge'
                roi=drawpolyline(ax);
                rmask=createMask(roi);
                Rfixed = imref2d(obj.Dimensions);
                tmask=imwarp(rmask,transform,'OutputView',Rfixed);
                Target=imgaussfilt(double(tmask),options.smoothing);

             case 'polyregion'
                roi=drawpolygon(ax);
                rmask=createMask(roi);
                Rfixed = imref2d(obj.Dimensions);
                tmask=imwarp(rmask,transform,'OutputView',Rfixed);
                Target=imgaussfilt(double(tmask),options.smoothing);
             otherwise
                disp('Invalid ROI type. Defaulting to points.')
                roi=drawpolyline(ax);
                i=roi.Position(:,1);
                j=roi.Position(:,2);
                [x,y]=transformPointsForward(transform,i,j);
                Target=obj.GenTargetFromSpots([y,x],obj.Dimensions,options.smoothing);
           end
           if options.append
                obj.Target = obj.Target + Target;
                obj.Target = obj.Target./max(obj.Target(:));
           else
                obj.Target = Target;
           end
        end
        
        function obj=SLM_Project(obj,options)
            arguments
                obj SLM_Device
                options.StackSelect double=0
            end
            if options.StackSelect==0
                Z_GS_loc=obj.Z_GS;
            else
                Z_GS_loc=obj.Hologram_Stack(options.StackSelect);
            end
               width = size(obj.Target,1);
               height = size(obj.Target,2);
               [X,Y] = meshgrid(1:width,1:height);
                
               %Creating a blazed diffraction grating in Y
               Z2 = angle(exp(-1i*(Y/obj.Blaze(1))))';

               %Creating a blazed diffraction grating in X
               Z1 = angle(exp(1i*(X/obj.Blaze(2))))';

               %Computes a phase for the focus
               focus = 0;%angle(exp(1i*((X-width/2).^2*obj.FocusX+(Y-height/2).^2*obj.FocusY)/2));

               %Sums all the previous elements
               SUM =  Z_GS_loc + Z1 + Z2 + focus;

               % Wraps to pi the contribution of all parts
               Z_SUM = wrapToPi(SUM);

               % Specifies the right range
               Z_TOTAL = 256*(Z_SUM+pi)/(2*pi);

               % Converts the computed phase into unit8 (as desired by the SLM)
               final_image = uint8(mod(Z_TOTAL,256));
            if obj.demo_mode==0
               calllib('Blink_C_wrapper', 'Write_image', final_image, 1);
            else
                figure
                imagesc(final_image)
                title('projected hologram')
            end
        end
    end

    methods(Static)

        function Target=GenTargetFromSpots(spotcoords,Dimensions,Width)
           [X,Y]=meshgrid(1:Dimensions(1),1:Dimensions(2));
           gauss = @(x0,y0,w0) exp(-((X-x0).^2+(Y-y0).^2)/w0.^2);
           Target=zeros(Dimensions);
           for i=1:size(spotcoords,1)
              Target=Target+gauss(spotcoords(i,1),spotcoords(i,2),Width)';
           end
        end

        function ClearSDK()
            if obj.demo_mode==0
               calllib('Blink_C_wrapper', 'Delete_SDK');
               unloadlibrary('Blink_C_wrapper');
            end
        end
    end

end
