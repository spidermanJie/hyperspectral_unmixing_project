%%   Generates  the  plots in Figure 7 for  SNR_SD 
%
%          (SNR_SD  == signal-to-noise-ratio spectral distribution)
%
%
%
%   Author: Jose Bioucas Dias (bioucas@lx.it.pt), November 2011
%
%%

   
clear all;
close all

%% load data set
load  SudP5SNR40
[B,n] = size(Y);

% estimate signal and noise correlation matrices
[w Rw] = estNoise(Y);
Rx = (Y-w)*(Y-w)'/n;

[U,S] = svd(Rx);
S = max(0,diag(S));

% compute SNR_SD
SNR_SD1 = S./diag(U'*Rw*U);


%% load data set
load  '..\..\datasets\SusgsP5SNR40'
[B,n] = size(Y);

% estimate signal and noise correlation matrices

[w Rw] = estNoise(Y);
Rx = (Y-w)*(Y-w)'/n;

[U,S] = svd(Rx);
S = max(0,diag(S));

SNR_SD2 = S./diag(U'*Rw*U) ;


%% load data set  Rcuprite

load  '..\..\datasets\Rcuprite'
[B,n] = size(Y);

[w Rw] = estNoise(Y);
Rx = (Y-w)*(Y-w)'/n ;

[U,S] = svd(Rx);
S = max(0,diag(S));

SNR_SD3 = S./diag(U'*Rw*U);

figure(1);
semilogy([SNR_SD1(1:50) SNR_SD2(1:50) SNR_SD3(1:50)], 'Linewidth',2)
title('SNR-SD')
xlabel('eigen direction')
legend('Sud','Susgs', 'Rcuprite' )
set(gca,'FontSize',16)
axis([0 50 1e-2 1e7])
set(gca,'YTick',[1e-2 1e0 1e2 1e4 1e6])

%print -depsc  snr_sd


