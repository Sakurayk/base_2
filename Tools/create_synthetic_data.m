%Author: Masaki Onuki (masaki.o@msp-lab.org)
function [Observed_mat,Miss_mat,A] = create_synthetic_data(h_size,w_size,Miss_rate,rank)

%A1 = ones(h_size/rank,w_size/rank)+0.5;

%A = blkdiag(A1,A1);
%for i=1:rank-2
%    A = blkdiag(A,A1);
%end
RGB = imread('Red-brick-wall-texture-3.jpg');
I = rgb2gray(RGB);
targetsize = [size(I,1),size(I,1)];
if rem(size(I,1),2)==1
    targetsize = [size(I,1)-1,size(I,1)-1];
end
rectangle = centerCropWindow2d(size(I),targetsize);
A = imcrop(I,rectangle);
A = double(I);
%A = abs(A-1); %ここを変えれば行けるのか?
w_size_miss = int32(w_size*Miss_rate);
Miss_pre_mat = [ones(h_size,w_size-w_size_miss),zeros(h_size,w_size_miss)];
Miss_vec = Miss_pre_mat(:);
Miss_mat = reshape(Miss_vec(randperm(h_size*w_size)),[h_size,w_size]);

Observed_mat = cast(double(A).*Miss_mat,'uint8');
