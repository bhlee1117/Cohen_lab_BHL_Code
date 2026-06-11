function [aligned_avgImg tform]=imReg_faster(avgImg,avgImg2)

% % Normalize input images
% fixed  = mat2gray(avgImg2);  % target image (reference)
% moving = mat2gray(avgImg);   % image to be aligned
% 
% % Automatically estimate rotation + translation
% tform = imregcorr(moving, fixed, 'rigid');
% 
% % Apply transformation to avgImg
% Rfixed = imref2d(size(fixed));
% aligned_avgImg = imwarp(moving, tform, 'OutputView', Rfixed);
% 
% % Display results
% nexttile([1 1]);
% imshowpair(aligned_avgImg, fixed, 'falsecolor');
% title('Aligned avgImg with avgImg2');

   fixed  = mat2gray(avgImg2);
    moving = mat2gray(avgImg);

    % Step 1: coarse alignment via phase correlation
    tform_init = imregcorr(moving, fixed, 'translation');

    % Step 2: refine with intensity-based registration
    optimizer = registration.optimizer.RegularStepGradientDescent;
    metric    = registration.metric.MeanSquares;
    optimizer.MaximumIterations = 300;
    optimizer.MinimumStepLength = 1e-5;

    tform = imregtform(moving, fixed, 'translation', optimizer, metric, ...
                       'InitialTransformation', tform_init);

    dx = tform.T(3,1);
    dy = tform.T(3,2);
    fprintf('Pixel shift to align: dx = %.2f px, dy = %.2f px\n', dx, dy);

    Rfixed = imref2d(size(fixed));
    aligned_avgImg = imwarp(moving, tform, 'OutputView', Rfixed);

    nexttile([1 1]);
    imshowpair(aligned_avgImg, fixed, 'falsecolor');
    title('Aligned');
end