clc
clearvars
close all
addpath Tools

%% 1:User Settings
imname = 'fallguys.jpg';
RGB = imread(imname);
%----The size of the data is determined----
%h_size = size(RGB,1);
%w_size = size(RGB,2);
w_size = 160;%The vertical size of the data (its default size is 1000)
h_size = 160;
target_size = [h_size,w_size];
%----The rank and missing rate for the data is determined----
rank = 10;% The rank of the data (10, 20, or  200 is used for the rank)
Miss_rate = 50;%The missing rate of the data (1, 10, or 20 is used for the missing rate)

Miss_percent = Miss_rate/100;

%----The synthetic data used in this experiment is created----
[V_RGB,L,Ground_Truth] = Synthesize_pic(imname,"./testing_mask_dataset/00017.png",target_size);

%----The used method is selected----
CPA_SVS = 1;% The CPA-based method is used if CPA_SVS = 1, or the exact method is used if CPA_SVS = 0.

if CPA_SVS == 1
    disp('The CPA-based method is selected for this experiment.')
elseif CPA_SVS == 0
    disp('The exact method is selected for this experiment.')
else
    Err_msg = sprintf('Error: You must select CPA_SVS = 1 or 0. Please, change it!');
    error(Err_msg);
    return;
end
disp(['default RMSE:',num2str(immse(V_RGB,Ground_Truth))])
%----If you select CPA_SVS = 1, you could change the approximation order for CPA----
Approx_order = 5;% Approx_order = 5, 10, 15, or 20 is used in this paper.
U_RGB = zeros(h_size,w_size,3);
time = zeros(1,3);
rmse = zeros(1,3);
[LoD,HiD] = wfilters('haar','d');
tstart = tic;

for k=1:3
    %% 2:Initialization of Some Variables
    V = V_RGB(:,:,k);
    U = U_RGB(:,:,k);
    Old_val = ones(h_size,w_size);
    z1 = ones(h_size,w_size);
    z2 = ones(h_size,w_size);
    z3 = ones(h_size,w_size);
    z4 = ones(h_size,w_size);
    d1 = ones(h_size,w_size);
    d2 = ones(h_size,w_size);
    d3 = ones(h_size,w_size);
    d4 = ones(h_size,w_size);

    stopcri = 1e-4; % stopping criterion
    maxiter = 1000; % maximum number of iteration

    I=ones(h_size,w_size);
    I_cheby=speye(w_size/2,w_size/2);

    max_level = 1;
    Q=zeros(w_size/2,h_size/2);
    Q_hat=zeros(w_size/2,h_size/2);
    %% 3:Optimization
    tic
    for i = 1:maxiter
        G=3*I+L;
        G=1./G;
        U = G.*((z1+z2+L.*z3+z4)-d1-d2-L.*d3-d4);


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
            h = @(x)(svd_kernel(x,5));
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
            D2=diag(sign(D2).*max(abs(D2)-0.5,0));
            z1 = S_l*D2*VT';
        end

        %prox:L1 norm
        B= U + d2;
        [cA,cH,cV,cD] = dwt2(B,LoD,HiD,'haar','spd');
        z2 = idwt2(cA,cH,cV,cD,'haar');

        %prox:Indicator
        z3=V;

        %prox:range
        C=U+d4;
        C(C<0)=0;
        C(C>1)=1;
        z4=C;

        %update
        d1 = d1+U-z1;
        d2 = d2+U -z2;
        d3 = d3+L.*U-z3;
        d4 = d4+ U - z4;

        Now_val = U;
        error(i) = sqrt(sum(sum((Now_val(:)-Old_val(:)).^2))) / sqrt(sum(sum((Now_val(:)).^2)));
        Old_val = Now_val;
        if error(i) < stopcri
            break;
        end
    end
    time(k) = toc;
    U_RGB(:,:,k) = U;
    rmse(k) = immse(U,Ground_Truth(:,:,k));
end
%% 4:Representation of the results
disp(['rmse:',num2str(immse(U_RGB,Ground_Truth))])
disp('--rmse--')
disp(rmse)
disp(['psnr:',num2str(psnr(U_RGB,Ground_Truth))])
disp('--time--')
disp(time)
disp(['total:',num2str(sum(time))])
disp('The results are shown.');
figure
subplot(2,2,1)
imshow(uint8(hsv2rgb(Ground_Truth)*255))
title({'Original data';})

subplot(2,2,2)
imshow(uint8(hsv2rgb(V_RGB)*255))
title({'Corrupted data';['[Corruption rate:', num2str(Miss_rate), '%]']})

if CPA_SVS == 1
    subplot(2,2,3)
    imshow(uint8(hsv2rgb(U_RGB)*255))
    title({'Resulting data (CPA-based method)';['[Approximation order:', num2str(Approx_order),']'];['[Number of iterations:', num2str(i),']']})
elseif CPA_SVS == 0
    subplot(2,2,3)
    imshow(uint8(hsv2rgb(U_RGB)*255))
    title({'Resulting data (Exact method)';['[Number of iterations:', num2str(i),']']})
end

subplot(2,2,4)
plot(1:length(error),error,'-','Linewidth',2);
title('Error')
xlabel('t')
ylabel('E_t')
ylim([-0.01 0.1])
grid on

    