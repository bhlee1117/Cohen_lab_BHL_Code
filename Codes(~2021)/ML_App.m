%% Scene Identification Using Image Data
% Image classification involves determining if an image contains some 
% specific object, feature, or activity. The goal of this example is to
% provide a strategy to construct a classifier that can automatically 
% detect which scene we are looking at
% This example uses function from the Computer Vision System Toolbox and
% Statistics and Machine Learning
% Copyright (c) 2017, MathWorks, Inc.

%% Description of the Data
% The dataset contains 4 scenes: Field, Auditorium, Beach ,Restaurant
% The images are photos of the scenes that have been taken from different
% angles, positions, and different lighting conditions. These variations make 
% this a challenging task.
%%%
%% 
% Data needs to be downloaded in order for this demo to run
% Please use your own data or download a database online
% Here's an example of a good source of scene data
% http://groups.csail.mit.edu/vision/SUN/

% Please note that all data in this example was resized to the same height
% and width dimensions. Montage will not work unless you resize all of the images
% This is done using the custom read function for ImageDatastore :
% readAndResizeImages.m
% You can alter the size of the images by opening that function

%% Load image data
% This assumes you have a directory: Scene_Categories
% with each scene in a subdirectory
clear
setDir  = 'E:\BACKUP\대학원\연구실\MY_Projects\In_vivo_imaging\Cell_detection\TXN_detection_deep_learning\Subtraction';
imds = imageDatastore(setDir,'IncludeSubfolders',true,'LabelSource',...
    'foldernames');       %^%      %#ok
%imds.ReadFcn = @readAndResizeImages;
%% Display Class Names and Counts
tbl = countEachLabel(imds)                                             %#ok
categories = tbl.Label;

%% Display Sampling of Image Data
visImds = splitEachLabel(imds,1,'randomize');

for ii = 1:2 % this assumes 4 categories of scenes
    subplot(2,2,ii);
    imagesc(visImds.readimage(ii));
    colormap('gray')
    title(char(visImds.Labels(ii)));
end

%% Pre-process Training Data: *Feature Extraction using Bag Of Words*
% Bag of features, also known as bag of visual words is one way to extract 
% features from images. To represent an image using this approach, an image 
% can be treated as a document and occurance of visual "words" in images
% are used to generate a histogram that represents an image.
%% Partition 700 images for training and 200 for testing
[training_set, test_set] = splitEachLabel(imds,0.7,'randomize');

%% Create Visual Vocabulary 
%extractorFcn = @exampleBagOfFeaturesExtractor;
tic 
bag = bagOfFeatures(imageSet(training_set.Files),'VocabularySize',500);%,'PointSelection','Grid'); 
scenedata = double(encode(bag,imageSet(training_set.Files))); 
toc 

%% Visualize Feature Vectors
imgs = training_set.splitEachLabel(1,'randomize',true);

img = readimage(imgs,1);
featureVector = encode(bag, img);
subplot(4,2,1); imshow(img);
subplot(4,2,2);
bar(featureVector);title('Visual Word Occurrences');xlabel('Visual Word Index');ylabel('Frequency');

img = readimage(imgs,2);
featureVector = encode(bag, img);
subplot(4,2,3); imshow(img);
subplot(4,2,4); 
bar(featureVector);title('Visual Word Occurrences');xlabel('Visual Word Index');ylabel('Frequency');
% 
% img = readimage(imgs,3);
% featureVector = encode(bag, img);
% subplot(4,2,5); imshow(img);
% subplot(4,2,6); 
% bar(featureVector);title('Visual Word Occurrences');xlabel('Visual Word Index');ylabel('Frequency');
% 
% img = readimage(imgs,4);
% featureVector = encode(bag, img);
% subplot(4,2,7); imshow(img);
% subplot(4,2,8); 
% bar(featureVector);title('Visual Word Occurrences');xlabel('Visual Word Index');ylabel('Frequency');

%% Create a Table using the encoded features
SceneImageData = array2table(scenedata);
%sceneType = categorical(repelem({training_set.Description}', [training_set.Count], 1));
SceneImageData.sceneType = training_set.Labels;

%% Create Visual Vocabulary 2
%extractorFcn = @exampleBagOfFeaturesExtractor;
tic 
extractorFcn=@customTXNExtractor;
bag = bagOfFeatures(imageSet(training_set.Files),...
    'CustomExtractor',extractorFcn,'StrongestFeatures',1,'VocabularySize',100);%,'PointSelection','Grid'); 
scenedata = double(encode(bag,imageSet(training_set.Files))); 
toc 

%% Visualize Feature Vectors
imgs = training_set.splitEachLabel(1,'randomize',true);

img = readimage(imgs,1);
featureVector = encode(bag, img);
subplot(4,2,1); imshow(img);
subplot(4,2,2);
bar(featureVector);title('Visual Word Occurrences');xlabel('Visual Word Index');ylabel('Frequency');

img = readimage(imgs,2);
featureVector = encode(bag, img);
subplot(4,2,3); imshow(img);
subplot(4,2,4); 
bar(featureVector);title('Visual Word Occurrences');xlabel('Visual Word Index');ylabel('Frequency');
% 
% img = readimage(imgs,3);
% featureVector = encode(bag, img);
% subplot(4,2,5); imshow(img);
% subplot(4,2,6); 
% bar(featureVector);title('Visual Word Occurrences');xlabel('Visual Word Index');ylabel('Frequency');
% 
% img = readimage(imgs,4);
% featureVector = encode(bag, img);
% subplot(4,2,7); imshow(img);
% subplot(4,2,8); 
% bar(featureVector);title('Visual Word Occurrences');xlabel('Visual Word Index');ylabel('Frequency');

%% Create a Table using the encoded features
SceneImageData = array2table(scenedata);
%sceneType = categorical(repelem({training_set.Description}', [training_set.Count], 1));
SceneImageData.sceneType = training_set.Labels;


%% Use the new features to train a model and assess its performance using 
classificationLearner

%% Test out accuracy on test set!
%montage(test_set.Files, 'Size', [1 7],'DisplayRange',[100 6000]);

testSceneData = double(encode(bag, test_set));
testSceneData = array2table(testSceneData,'VariableNames',trainedModel.RequiredVariables);

actualSceneType = test_set.Labels;

predictedOutcome = trainedModel.predictFcn(testSceneData);
correctPredictions = (predictedOutcome == actualSceneType);
validationAccuracy = sum(correctPredictions)/length(predictedOutcome) %#ok

%% Visualize how the classifier works
% figure(2);
% random_num = randi(length(test_set.Labels));
% 
% img = test_set.readimage(random_num);
% 
% imshow(img)
% % Add code here to invoke the trained classifier
% imagefeatures = double(encode(bag, img));
% % Find two closest matches for each feature
% [bestGuess, score] = predict(trainedClassifier.ClassificationSVM,imagefeatures);
% % Display the string label for img
% if bestGuess == test_set.Labels(random_num)
% 	titleColor = [0 0.8 0];
% else
% 	titleColor = 'r';
% end
% title(sprintf('Best Guess: %s; Actual: %s',...
% 	char(bestGuess),test_set.Labels(random_num)),...
% 	'color',titleColor)