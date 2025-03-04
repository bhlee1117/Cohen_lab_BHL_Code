%% 19-06-2020 :  This function was written in the case the im_straighter 
%% algorithm is not working properly. The following algoritm is also calculating
%% the in-focus image but after defining manually the planes. 
%% ==========================================================

clear all
close all
clc

WindowSize = 512;
Current_analyzed_folder = uigetdir('C:\Users\sCMOS-1\Desktop_');

%% Ask the number of ROIs that were acquired during each cycle
%% -----------------------------------------------------------

prompt = {'Enter the number of ROIs:'};
title = 'Input';
dims = [1 35];
definput = {'9'};
answer = inputdlg(prompt,title,dims,definput);

N_ROI = str2double(answer{1});

%% Look for the tiff files and start the calculation of the in-focus image
%% =======================================================================

[TIFF_FileName, ~] = Look_For_TIFF_Files_dcimg_conversion(Current_analyzed_folder);
N_tiff = size(TIFF_FileName,1);

for n_roi = 1 : N_ROI
    
    Selected_files = zeros(N_tiff,1);
    Char = strcat('ROI_', num2str(n_roi));
    
    % Select the files according to the ROI number
    % --------------------------------------------
    
    for n_file = 1 : N_tiff
        if ~isempty(strfind(TIFF_FileName{n_file}, Char))
            Selected_files(n_file) = 1;
        end
    end
    
    TIFF_selected = TIFF_FileName(Selected_files==1);
    
    if ~isempty(TIFF_selected)
        In_Focus_dir = strcat(Current_analyzed_folder, '\In_Focus_images');
        mkdir(In_Focus_dir)
        In_Focus_saving_folder = strcat(In_Focus_dir, '\ROI_', num2str(n_roi));
        mkdir(In_Focus_saving_folder)
        
        if n_roi == 1
            
            % Define the parameters for the calculation of the planes
            % -------------------------------------------------------
            
            cd(Current_analyzed_folder)
            imName = TIFF_selected{1};
            ImInfo = imfinfo(imName);
            Ly = ImInfo(1).Width;
            Lx = ImInfo(1).Height;
            
            if Lx == Ly
                NROItot = Lx/WindowSize;
            else
                warndlg('The images are expected to be square. The calculation is aborted')
                uiwait(warndlg)
                delete(warndlg)
                return
            end
            
            NPlanes = uint8(size(ImInfo,1));
            
            % Define manually the reference planes for the first image
            % --------------------------------------------------------
            
            Ref_plane = Define_ref_planes(imName, WindowSize, NROItot, NPlanes);
        end
        
        % Calcule the in_focus image according to the the calculated
        % reference plane and save the new image
        % --------------------------------------
        
        template = zeros(Lx,Ly);
        N_TIFF = size(TIFF_selected,1);
        
        for n_file = 1 : N_TIFF
            
            imName = TIFF_selected{n_file};
            cd(Current_analyzed_folder)
            
            in_focus_im = Apply_ref_planes(imName, Ref_plane, WindowSize, NROItot, template);
            
            % The new "in-focus" image is then saved in a folder with a name
            % --------------------------------------------------------------
            
            cd(In_Focus_saving_folder)
            t = Tiff(imName, 'w');
            
            tagstruct = struct('ImageLength', size(in_focus_im,1), ...
                'ImageWidth', size(in_focus_im,2), ...
                'BitsPerSample', 16, ...
                'Photometric', Tiff.Photometric.MinIsBlack, ...
                'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
            
            t.setTag(tagstruct);
            t.write(uint16(in_focus_im));
            t.close();
            
            disp(strcat('Image #', num2str(n_file,'%03d'), ' is saved'))
        end
    end
end

%% Manual selection of the in-focus plane
%% =======================================

function Ref_plane = Define_ref_planes(imName, WindowSize, NROItot, NPlanes)

% AllPlanes = zeros(size(tiff,1),NROItot^2);
% Nelement = (WindowSize*WindowSize)-1;

figure('units','normalized','outerposition',[0 0 1 1])
N_plot = ceil(sqrt(double(NPlanes)));
N_plot = double(N_plot);
Ref_plane = zeros(NROItot,NROItot);

for ROI_x = 1:1:NROItot
    for ROI_y = 1:1:NROItot
        
        Rect = [(ROI_x-1)*WindowSize+1,(ROI_y-1)*WindowSize+1,WindowSize-1,WindowSize-1];
        for plane = 1 : NPlanes
            im(:,:,plane) = imread(imName,plane);
            ImCrop = imcrop(im(:,:,plane), Rect);
            
            subplot(N_plot,N_plot,double(plane))
            imagesc(ImCrop)
            axis image
            axis off
            colormap('Gray')
            title(strcat('plane #', num2str(plane)))
        end
        
        prompt = {'Enter the in-focus plane:'};
        dlg_title = 'Input';
        num_lines = 1;
        defaultans = {''};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
        
        Ref_plane(ROI_x,ROI_y) = str2double(answer{1});
    end
end

figure
imagesc(Ref_plane)
axis image
colorbar

end

%% Calculate the in_focus image
%% =============================

function template = Apply_ref_planes(imName, Ref_plane, WindowSize, NROItot, template)

Planes = unique(Ref_plane);

for nplane = 1 : size(Planes,1)
    
    Selected_plane = Planes(nplane);
    im = imread(imName,Selected_plane);
    
    for row = 1 : NROItot
        for col = 1 : NROItot
            if Ref_plane(row, col) == Selected_plane
                Rect = [(row-1)*WindowSize+1,(col-1)*WindowSize+1,WindowSize-1,WindowSize-1];
                ImCrop = imcrop(im, Rect);
                template((col-1)*WindowSize+1:(col-1)*WindowSize+WindowSize , (row-1)*WindowSize+1:(row-1)*WindowSize+WindowSize) = ImCrop;
            end
        end
    end
end
end
