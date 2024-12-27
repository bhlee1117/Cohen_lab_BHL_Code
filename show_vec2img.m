function show_vec2img(ImgStack,V)

% Sample Data
N = size(ImgStack,3); % Number of images
X = size(ImgStack,1); % Image height
Y = size(ImgStack,2); % Image width

% Normalize the vector to [0, 1] for coloring (if not already normalized)
V_norm = (V - min(V)) / (max(V) - min(V))+0.1;

ProjectedImage = zeros(X, Y);

% Combine slices by weighting with vector values
ProjectedImage = max(ImgStack.*reshape(V_norm,1,1,[]),[],3);
imagesc(ProjectedImage)
end