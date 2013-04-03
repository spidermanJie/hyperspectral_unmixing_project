%%    
%  Projects spectra on the signal subspace and plots Fig. 9 
%
%
%   Author: Jose Bioucas Dias (bioucas@lx.it.pt),  November 2011
%
%%   
   
clear all;
close all

verbose = 'on';

%% load data set  SusgsP5SNR30
load  SusgsP5SNR30
[B,n] = size(Y);

% estimate  noise
[w Rw] = estNoise(Y, verbose);
% estimate subspace
[kf,Ek] = hysime(Y,w,Rw,verbose);

% estimate unprojected SNR
SNR_sim_unp = trace(Y*Y'/n - Rw)/trace(Rw);
SNR_sim_unp = 10*log10(SNR_sim_unp);

% estimate projected SNR
SNR_sim_proj = trace(Y*Y'/n - Rw)/trace(Ek'*Rw*Ek);
SNR_sim_proj = 10*log10(SNR_sim_proj);


figure(1)
plot([M*x(:,1000) Y(:,1000) Ek*Ek'*Y(:,1000)],'LineWidth',2)
title('Simulated data set','FontSize',16)
xlabel('band','FontSize',16)
legend('original','noisy (SNR = 30 dB)', 'projected (SNR = 46.6 dB)')
set(gca,'FontSize',16)

%print -depsc  projected_simulated


%% load data set  Rcuprite

load  'Rcuprite'
[B,n] = size(Y);

% estimate noise
[w Rw] = estNoise(Y);


% estimate subspace
[kf,Ek] = hysime(Y,w,Rw,verbose);

% estimate unprojected SNR
SNR_cuprite_unp = trace(Y*Y'/n - Rw)/trace(Rw);
SNR_cuprite_unp = 10*log10(SNR_cuprite_unp);

% estimate projected SNR
SNR_cuprite_proj = trace(Y*Y'/n - Rw)/trace(Ek'*Rw*Ek);
SNR_cuprite_proj = 10*log10(SNR_cuprite_proj);

figure(2)
plot([Y(:,1000) Ek*Ek'*Y(:,1000)],'LineWidth',2)
title('Cuprite data set','FontSize',16)
xlabel('band','FontSize',16)
legend('noisy (SNR = 42.5 dB)', 'projected (SNR = 47.5 dB)' )
set(gca,'FontSize',16)

%print -depsc  projected_cuprite



