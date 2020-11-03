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
%----The size of the data is determined----
w_size = 128;%The vertical size of the data (its default size is 1000)
h_size = 128;
RGB = imread('modern-hexagonal-glowing-blue-medical-background-texture-pattern-honeycombs-different-level-d-rendering-illustration-futuristic-165624902.jpg');
I_2 = rgb2gray(RGB);
[I_R,I_G,I_B] = imsplit(RGB);
%----The size of the data is determined----
%h_size = size(I_2,1);
%w_size = size(I_2,2);
%----The rank and missing rate for the data is determined----
rank = 10;% The rank of the data (10, 20, or  200 is used for the rank)
Miss_rate = 30;%The missing rate of the data (1, 10, or 20 is used for the missing rate)

Miss_percent = Miss_rate/100;

%----The synthetic data used in this experiment is created----
[V_RGB,L,Ground_Truth]=create_synthetic_data(h_size,w_size,Miss_percent,rank);

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
U_RGB = zeros(h_size,w_size,3);

for k=1:3
    %% 2:Initialization of Some Variables
    V = V_RGB(:,:,k);
    U = U_RGB(:,:,k);
    Old_val = ones(h_size,w_size);
    z1 = ones(h_size,w_size);
    z2 = dct2(ones(h_size,w_size));
    z3 = ones(h_size,w_size);
    z4 = ones(h_size,w_size);
    d1 = ones(h_size,w_size);
    d2 = dct2(ones(h_size,w_size));
    d3 = ones(h_size,w_size);
    d4 = ones(h_size,w_size);

    stopcri = 1e-4; % stopping criterion
    maxiter = 200; % maximum number of iteration

    I=ones(h_size,w_size);
    I_cheby=speye(w_size/2,w_size/2);

    max_level = 1;
    Q=zeros(w_size/2,h_size/2);
    Q_hat=zeros(w_size/2,h_size/2);

    %% 3:Optimization
    disp('------------------------------Optimization------------------------------')
    disp('ADMM is running...');

    for i = 1:maxiter
        G=3*I+L;
        G=1./G;
        U = G.*((z1+idct2(z2)+L.*z3+z4)-d1-idct2(d2)-L.*d3-d4);


        %prox:nuclear norm
        if CPA_SVS == 1
            %----The singular value shrinkage using the CPA-based method----
            A=U+d1;
            A_1 = A'*A;
            [Q,~,~,~]=dwt2(A_1,'haar');
            d=eigs(Q,1);
            if d==0
                d=100000;
            end
            h = @(x)(svd_kernel(x,6));
            Range = [0 d];
            Coeff=chebyshev_coefficient(h, Approx_order,Range);%Chebyshev coefficients are derived
            Q_hat = chebyshev_oprator(I_cheby,Q,Coeff,Range);%Singular value shrinkage is performed by using CPA

            [A_hat]=idwt2(Q_hat,zeros(size(Q_hat)),zeros(size(Q_hat)),zeros(size(Q_hat)),'haar');
            A_hat = A*A_hat;
            z1 = A_hat;
        elseif CPA_SVS == 0
            %----The singular value shrinkage using the exact method----
            A=U+d1;
            [S_l,D,VT]=svd(A);
            D2=diag(D);
            D2=diag(sign(D2).*max(abs(D2)-6,0));
            z1 = S_l*D2*VT';
        end

        %prox:L1 norm
        B=dct2(U)+d2;
        z2 = sign(B).*max(abs(B)-0.1,0);

        %prox:Indicator
        z3=V;

        %prox:range
        C=U+d4;
        C(C<0)=0;
        C(C>1)=1;
        z4=C;

        %update
        d1 = d1+U-z1;
        d2 = d2+dct2(U)-z2;
        d3 = d3+L.*U-z3;
        d4 = d4+ U - z4;

        Now_val = U;
        error(i) = sqrt(sum(sum((Now_val(:)-Old_val(:)).^2))) / sqrt(sum(sum((Now_val(:)).^2)));
        Old_val = Now_val;
        if error(i) < stopcri
            break;
        end
    end
    U_RGB(:,:,k) = U;
end
%% 4:Representation of the results
disp('------------------------------Results------------------------------')
disp('The results are shown.');
figure
subplot(2,2,1)
imshow(cast(Ground_Truth,'uint8'))
title({'Original data';['[Matrix rank:', num2str(rank),']']})

subplot(2,2,2)
imshow(cast(V_RGB,'uint8'))
title({'Corrupted data';['[Corruption rate:', num2str(Miss_rate), '%]']})

if CPA_SVS == 1
    subplot(2,2,3)
    imshow(cast(U_RGB,'uint8'))
    title({'Resulting data (CPA-based method)';['[Approximation order:', num2str(Approx_order),']'];['[Number of iterations:', num2str(i),']']})
elseif CPA_SVS == 0
    subplot(2,2,3)
    imshow(cast(U_RGB,'uint8'))
    title({'Resulting data (Exact method)';['[Number of iterations:', num2str(i),']']})
end

subplot(2,2,4)
plot(1:length(error),error,'-','Linewidth',2);
title('Error')
xlabel('t')
ylabel('E_t')
ylim([-0.01 0.1])
grid on

    