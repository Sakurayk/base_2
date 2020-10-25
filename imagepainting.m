%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo: Fast singular value shrinkage using the exact method for synthetic data
% Note: This exmeriments were conducted at Section V-F in the following
%           paper.
%
% Author: Masaki Onuki (masaki.o@msp-lab.org)
% Last version: Aug 17, 2017
% Article: M. Onuki, S. Ono, K. Shirai, Y. Tanaka,
% "Fast Singular Value Shrinkage with Chebyshev Polynomial Approximation Based on Signal Sparsity,"
% IEEE Transactions on Signal Processing (submitted).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clearvars
close all
addpath Tools

%% 1:User Settings
%image reading
RGB = imread('Red-brick-wall-texture-3.jpg');
I = rgb2gray(RGB);

%square image
targetsize = [size(I,1),size(I,1)];
if rem(size(I,1),2)==1
    targetsize = [size(I,1)-1,size(I,1)-1];
end
rectangle = centerCropWindow2d(size(I),targetsize);
I_2 = imcrop(I,rectangle);
I_2 = double(I_2);

%----The size of the data is determined----
h_size = size(I_2,1);
w_size = size(I_2,2);
%h_size = 1000;%The horizontal size of the data (its default size is 1000)
%w_size = 1000;%The vertical size of the data (its default size is 1000)

%----The rank and missing rate for the data is determined----
rank = 10;% The rank of the data (10, 20, or  200 is used for the rank)
Miss_rate = 20;%The missing rate of the data (1, 10, or 20 is used for the missing rate)

Miss_percent = Miss_rate/100;

%----The synthetic data used in this experiment is created----
[Observed_image,Miss_mat_L,Original_data_I]=create_synthetic_data_2(h_size,w_size,Miss_percent,rank,I_2);

%----The used method is selected----
CPA_SVS = 0;% The CPA-based method is used if CPA_SVS = 1, or the exact method is used if CPA_SVS = 0.

if CPA_SVS == 1
    disp('The CPA-based method is selected for this experiment.')
elseif CPA_SVS == 0
    disp('The exact method is selected for this experiment.')
else
    Err_msg = sprintf('Error: You must select CPA_SVS = 1 or 0. Please, change it!');
    error(Err_msg);
    return;
end

%----If you select CPA_SVS = 1, you could change the approximation order for CPA----
Approx_order = 20;% Approx_order = 5, 10, 15, or 20 is used in this paper.


%% 2:Initialization of Some Variables
omega = Miss_mat_L;
omega_bar = Miss_mat_L==0;
M_delta = zeros(h_size,w_size);
M_delta(495:755,495:755) = 1;
M_delta = M_delta - omega_bar;

%あとでベクトル化する

I_vec = Original_data_I;
L_vec = Miss_mat_L;
M_vec = Observed_image;
dct_observed = dct2(Observed_image);

z1 = eye(h_size,w_size);
z2 = dct_observed;
z3 = Miss_mat_L;%omega
z4 = eye(h_size,w_size);
z5 = Miss_mat_L ==0;%omega_bar
K = [z1;z2;z3;z4;z5];
z = [z1*Observed_image;z2*Observed_image;z3*Observed_image;z4*Observed_image;z5*Observed_image];

u1 = zeros(h_size,w_size);
u2 = zeros(h_size,w_size);
u3 = zeros(h_size,w_size);
u4 = zeros(h_size,w_size);
u5 = zeros(h_size,w_size);
u = [u1',u2',u3',u4',u5']';

stopcri = 1e-4; % stopping criterion
maxiter = 20; % maximum number of iteration

I=ones(h_size,w_size);
I_cheby=speye(int32(w_size/2),int32(h_size/2));

max_level = 1;
Q=zeros(int32(w_size/2),int32(h_size/2));
Q_hat=zeros(int32(w_size/2),int32(h_size/2));
Old_val = u;

%% 3:Optimization
disp('------------------------------Optimization------------------------------')
disp('ADMM is running...');

for i = 1:maxiter
    %G=3*I+Miss_mat_L;
    %G=1./G;
    L_vec_ans = (K.'*K)\K.'*(z-u);
    %prox:nuclear norm
    if CPA_SVS == 1
        %----The singular value shrinkage using the CPA-based method----
        %A=U+d1;
        %A_1 = A'*A;
        A_z_1 = L_vec_ans + u1;
        A_1 = A_z_1'*A_z_1;
        [Q,~,~,~]=dwt2(A_1,'haar');
        d=eigs(Q,1);
        if d==0
            d=100000;
        end
        h = @(x)(svd_kernel(x,5));
        Range = [0 d];
        Coeff=chebyshev_coefficient(h, Approx_order,Range);%Chebyshev coefficients are derived
        Q_hat = chebyshev_oprator(I_cheby,Q,Coeff,Range);%Singular value shrinkage is performed by using CPA
        
        [A_hat]=idwt2(Q_hat,zeros(size(Q_hat)),zeros(size(Q_hat)),zeros(size(Q_hat)),'haar');
        A_hat = A_1*A_hat;
        z1 = A_hat;
                
    elseif CPA_SVS == 0
        %----The singular value shrinkage using the exact method----
        A_z_1 = L_vec_ans + u1;
        [S_l,D,VT]=svd(A_z_1);
        D2=diag(D);
        D2=diag(sign(D2).*max(abs(D2)-0.1,0));
        z1 = S_l*D2*VT';
    end
    
    %l1norm  
    B=dct_observed*L_vec_ans+u2;
    z2 = sign(B).*max(abs(B)-0.1,0);

    
    %prox:Indicator
    z3 = Observed_image;
    %元のピクセルを維持する
    
    %prox:range
    C = L_vec_ans + u4;
    C(C<=0)=0;
    C(C>=1)=1;
    z4 = C;
    
    M_val = mean(omega_bar*L_vec_ans+u5,'all')-sum((M_delta*Observed_image),'all')/(sum(M_delta,'all'));
    z5 = z5+omega_bar*M_val;
    
    u = u + K*L_vec_ans - [z1;z2;z3;z4;z5];
    u1 = u(1:h_size,:);
    u2 = u(h_size+1:2*h_size,:);
    u3 = u(2*h_size+1:3*h_size,:);
    u4 = u(3*h_size+1:4*h_size,:);
    u5 = u(4*h_size+1:5*h_size,:);
    %update
    
    Now_val = u ;
    error(i) = sqrt(sum(sum((Now_val(:)-Old_val(:)).^2))) / sqrt(sum(sum((Now_val(:)).^2)));
    Old_val = Now_val;
    if error(i) < stopcri
        break;
    end
end

%% 4:Representation of the results
disp('------------------------------Results------------------------------')
disp('The results are shown.');
figure
subplot(2,2,1)
imshow(uint8(Original_data_I))
title({'Original data';['[Matrix rank:', num2str(rank),']']})

subplot(2,2,2)
imshow(uint8(Observed_image))
title({'Corrupted data';['[Corruption rate:', num2str(Miss_rate), '%]']})

if CPA_SVS == 1
    subplot(2,2,3)
    result_image = reshape(L_vec_ans,[h_size,w_size]);
    imshow(uint8(L_vec_ans))
    title({'Resulting data (CPA-based method)';['[Approximation order:', num2str(Approx_order),']'];['[Number of iterations:', num2str(i),']']})
elseif CPA_SVS == 0
    subplot(2,2,3)
    imshow(uint8(L_vec_ans))
    title({'Resulting data (Exact method)';['[Number of iterations:', num2str(i),']']})
end

subplot(2,2,4)
plot(1:length(error),error,'-','Linewidth',2);
title('Error')
xlabel('t')
ylabel('E_t')
ylim([-0.01 1.0])
grid on

    