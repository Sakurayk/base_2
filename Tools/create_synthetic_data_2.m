function [Observed_mat,Miss_mat,A] = create_synthetic_data_2(h_size,w_size,Miss_rate,rank,read_image)

%A1 = ones(h_size/rank,w_size/rank)+0.5;
%A = blkdiag(A1,A1);
%for i=1:rank-2
%    A = blkdiag(A,A1);
%end

%A = abs(A-1); %ここを変えれば行けるのか?
A = read_image;
%RGB2HSV
A = rgb2hsv(A);
A = double(A);
[A_1,A_2,A_3] = imsplit(A);
%Miss_pre_mat = [ones(h_size,w_size-w_size_miss),zeros(h_size,w_size_miss)];
%Miss_vec = Miss_pre_mat(:);
%Miss_mat = reshape(Miss_vec(randperm(h_size*w_size)),[h_size,w_size]);
Miss_mat = ones(h_size,w_size);
Miss_mat(25:35,25:35) = 0;

Observed_mat_1 = double(A_1).*Miss_mat;
Observed_mat_2 = double(A_2).*Miss_mat;
Observed_mat_3 = double(A_3).*Miss_mat;
Observed_mat = cat(3,Observed_mat_1,Observed_mat_2,Observed_mat_3);
