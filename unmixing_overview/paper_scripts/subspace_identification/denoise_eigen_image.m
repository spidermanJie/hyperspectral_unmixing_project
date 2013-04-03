%% denoise a low SNR eigen-image
%
%   Uses the BM3D  algorithm  
%
%    K. Dabov, A. Foi, V. Katkovnik, and K. Egiazarian, “Image denoising by
%    sparse 3D transform-domain collaborative filtering,” IEEE Trans. 
%    Image Process., vol. 16, no. 8, pp. 2080-2095, August 2007.
%
%   to denoise noisy eigen-images
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT NOTE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      This script uses the package BM3D  (v1.9 (26 August 2011))
%      to implement the denoising algorithm BM3D introduced in 
%
%      K. Dabov, A. Foi, V. Katkovnik, and K. Egiazarian, "Image denoising by
%      sparse 3D transform-domain collaborative filtering," IEEE Trans. 
%      Image Process., vol. 16, no. 8, pp. 2080-2095, August 2007.
%
%      The BM3D package  is available at the 
%      BM3D web page:  http://www.cs.tut.fi/~foi/GCF-BM3D
%
%       Download this package and install it is the folder include/BM3D
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Jose Bioucas Dias (bioucas@lx.ir.pt), Nobember, 2011
%
%%
clear all;
close all

%% load data set  Rcuprite

load  Rcuprite
[B,n] = size(Y);

[w Rw] = estNoise(Y);
Rx = (Y-w)*(Y-w)'/n ;

[E, S] = svd(Rx); 


% choose band

for i=[1:18]
    % produce eigen-image
    eigen_im = E(:,i)'*Y;
    min_x = min(eigen_im);    
    max_x = max(eigen_im);
    eigen_im = eigen_im - min_x;
    scale = max_x-min_x;
    
    %scale to [0,1]
    eigen_im = reshape(eigen_im, Lines, Columns)/scale;
    figure(1);
    imagesc(eigen_im);
    
    %estimate noise from Rw
    sigma = sqrt(E(:,i)'*Rw*E(:,i))/scale;
    
    % denoise  with BM3D
    [dummy, filt_eigen_im] = BM3D(1,eigen_im, sigma*255);

    figure(2);
    imagesc(filt_eigen_im);
    
    % show the the diffrence between the filtered and the noise
    % eigendirection
    figure(3);
    imagesc(filt_eigen_im - eigen_im);
    colorbar;
    drawnow;
    
    Yf(i,:) = reshape(filt_eigen_im*scale + min_x, 1,Lines*Columns);
  
end


%% plot images


figure(1);
imagesc(eigen_im*scale + min_x);
axis off
title('Noisy eigen-image no. 18', 'FontSize', 16)
colorbar
   
figure(2);
imagesc(filt_eigen_im*scale + min_x);
axis off
title('Filtered eigen-image no. 18', 'FontSize', 16)
colorbar

figure(3);
imagesc((filt_eigen_im - eigen_im)*scale );
axis off
title('(Noisy - Filtered) eigen-image 18', 'FontSize', 16)
colorbar

index = 17;

% scatter plots (noisy and filterd) 
figure(4)
plot( E(:,index)'*Y(:,1:10:end), E(:,index+1)'*Y(:,1:10:end), '.',Yf(index,1:10:end), Yf(index+1,1:10:end), '.')
title('scatter plots of eigen-images','FontSize', 16)
xlabel('eigen-image 17','FontSize', 16)
ylabel('eigen-image 18','FontSize', 16)
axis([-0.05, 0.1, -0.03 0.04]);
legend('noisy', 'denoised')
set(gca, 'FontSize', 16)





