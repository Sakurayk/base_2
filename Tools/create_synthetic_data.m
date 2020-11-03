%Author: Masaki Onuki (masaki.o@msp-lab.org)
function [Observed_mat,Miss_mat,A] = create_synthetic_data(h_size,w_size,Miss_rate,rank)

%A1 = ones(h_size/rank,w_size/rank)+0.5;

%A = blkdiag(A1,A1);
%for i=1:rank-2
%    A = blkdiag(A,A1);
%end
RGB = imread('modern-hexagonal-glowing-blue-medical-background-texture-pattern-honeycombs-different-level-d-rendering-illustration-futuristic-165624902.jpg');
targetsize = [h_size,w_size];
rectangle = centerCropWindow2d([size(RGB,1),size(RGB,2)],targetsize);
RGB = imcrop(RGB,rectangle);
A = double(RGB);
%RGB2HSV
%A = rgb2hsv(RGB);
%A = double(A);
[A_R,A_G,A_B] = imsplit(A);
%A = abs(A-1); %ここを変えれば行けるのか?
w_size_miss = int32(w_size*Miss_rate);
Miss_pre_mat = [ones(h_size,w_size-w_size_miss),zeros(h_size,w_size_miss)];
Miss_vec = Miss_pre_mat(:);
Miss_mat = reshape(Miss_vec(randperm(h_size*w_size)),[h_size,w_size]);

Observed_mat_R = double(A_R).*Miss_mat;
Observed_mat_G = double(A_G).*Miss_mat;
Observed_mat_B = double(A_B).*Miss_mat;
Observed_mat = cat(3,Observed_mat_R,Observed_mat_G,Observed_mat_B);


