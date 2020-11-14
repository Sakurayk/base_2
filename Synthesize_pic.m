function [outputfile,Miss_mat,HSV] = Synthesize_pic(base_image_name,add_pic_name,targetsize)
% synthesize image for image inpainting

%read image
I = imread(base_image_name);
rectangle = centerCropWindow2d([size(I,1),size(I,2)],targetsize);
RGB = imcrop(I,rectangle);
HSV = rgb2hsv(RGB);
HSV = double(HSV);
[I_R,I_G,I_B] = imsplit(HSV);

%read synthesize image
I_2 = imread(add_pic_name);
Miss_mat = imbinarize(imresize(I_2,targetsize));
Miss_mat = (Miss_mat == 0); %binary inverse

Observed_mat_R = double(I_R).*Miss_mat;
Observed_mat_G = double(I_G).*Miss_mat;
Observed_mat_B = double(I_B).*Miss_mat;
outputfile = cat(3,Observed_mat_R,Observed_mat_G,Observed_mat_B);

end

