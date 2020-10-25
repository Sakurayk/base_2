%Author: Masaki Onuki (masaki.o@msp-lab.org)
function [Observed_mat,Miss_mat,A] = create_synthetic_data_2(h_size,w_size,Miss_rate,rank,read_image)

%A1 = ones(h_size/rank,w_size/rank)+0.5;

%A = blkdiag(A1,A1);
%for i=1:rank-2
%    A = blkdiag(A,A1);
%end

%A = abs(A-1); %ここを変えれば行けるのか?
A = read_image;
w_size_miss = int32(w_size*Miss_rate);
%Miss_pre_mat = [ones(h_size,w_size-w_size_miss),zeros(h_size,w_size_miss)];
%Miss_vec = Miss_pre_mat(:);
%Miss_mat = reshape(Miss_vec(randperm(h_size*w_size)),[h_size,w_size]);
Miss_mat = ones(h_size,w_size);
Miss_mat(500:750,500:750) = 0;

Observed_mat = double(A).*Miss_mat;
