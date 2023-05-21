function im = loadImage (filename)
image_name = dir([filename, '.tif']);
if ~isempty(image_name)
    t = Tiff([image_name.folder, '/',image_name.name],'r');
    im = read(t);
else
    image_name = dir([filename, '.mat']);
    load (image_name.name)
    im = andor_data.imageData;
end