clc
clearvars
close all
addpath Tools

%% User Settings
%----The size of the data is determined----
h_size = 160;
w_size = 160;%The vertical size of the data (its default size is 1000)
CPA_SVS = 1;
RGB = imread('modern-hexagonal-glowing-blue-medical-background-texture-pattern-honeycombs-different-level-d-rendering-illustration-futuristic-165624902.jpg');
I_2 = rgb2gray(RGB);
%----The size of the data is determined----
%h_size = size(I_2,1);
%w_size = size(I_2,2);
rank = 10;% The rank of the data (10, 20, or  200 is used for the rank)

%----The rank and missing rate for the data is determined----
Miss_rate = 10;%The missing rate of the data (1, 10, or 20 is used for the missing rate)

Miss_percent = Miss_rate/100;
%----The synthetic data used in this experiment is created----
[V,L,Ground_Truth]=create_synthetic_data(h_size,w_size,Miss_percent,rank);

%----split RGB
[V_R,V_G,V_B] = imsplit(V);

%% optimization
[U_R,R_error] = optimization(h_size,w_size,V_R,L);
[U_G,G_error] = optimization(h_size,w_size,V_G,L);
[U_B,B_error] = optimization(h_size,w_size,V_B,L);

%----concat
U = cat(3,U_R,U_G,U_B);
error_max = max([length(R_error),length(G_error),length(B_error)]);

%% Representation of the results
disp('------------------------------Results------------------------------')
disp('The results are shown.');
figure
subplot(2,2,1)
imshow(cast(RGB,'uint8'))
title({'Original data';})

subplot(2,2,2)
imshow(cast(hsv2rgb(V),'uint8'))
title({'Corrupted data';['[Corruption rate:', num2str(Miss_rate), '%]']})

if CPA_SVS == 1
    subplot(2,2,3)
    imshow(cast(hsv2rgb(U),'uint8'))
    title({'Resulting data (CPA-based method)';})
elseif CPA_SVS == 0
    subplot(2,2,3)
    imshow(U)
    title({'Resulting data (Exact method)';['[Number of iterations:', num2str(i),']']})
end

subplot(2,2,4)
plot(1:error_max,R_error,'-.r','Linewidth',2);
plot(1:error_max,G_error,'-.g','Linewidth',2);
plot(1:error_max,B_error,'-.b','Linewidth',2);
title('Error')
xlabel('t')
ylabel('E_t')
ylim([-0.01 0.1])
grid on

    