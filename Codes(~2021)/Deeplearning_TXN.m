%% Train, Evaluate, and Apply Image Category Classifier
%
clear
%% 
% Load two image categories.
setDir  = 'E:\BACKUP\대학원\연구실\MY_Projects\In_vivo_imaging\Cell_detection\TXN_detection_deep_learning\Subtraction';
imds = imageDatastore(setDir,'IncludeSubfolders',true,'LabelSource',...
    'foldernames');
%% 
% Split the data set into a training and test data. Pick 30% of images 
% from each set for the training data and the remainder  70% for the 
% test data.
[trainingSet,testSet] = splitEachLabel(imds,0.85,'randomize');
%% 
% Create bag of visual words.
bag = bagOfFeatures(trainingSet);
%% 
classificationLearner

%%
% Train a classifier with the training sets.
categoryClassifier = trainImageCategoryClassifier(trainingSet,bag);
%% 
% Evaluate the classifier using test images. Display the confusion matrix.
confMatrix = evaluate(categoryClassifier,testSet)
%% 
% Find the average accuracy of the classification.
mean(diag(confMatrix))
%% 
% Apply the newly trained classifier to categorize new images.

for i=1:7
img = imread(fullfile('E:\BACKUP\대학원\연구실\MY_Projects\In vivo imaging\Cell_detection\TXN_detection_deep_learning\Sub_Test\Test\',['Sub_' num2str(i) '.tif.tif']));
[labelIdx, score] = predict(categoryClassifier,img);
result(i,1)=[labelIdx];
end
%% 
% Display the classification label.
categoryClassifier.Labels(labelIdx)