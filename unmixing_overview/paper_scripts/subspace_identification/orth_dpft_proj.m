%%   Computes the orthogonal the DPFT projections and the angle with  the 
%    unprojected data
%
%
%
%   Author: Jose Bioucas Dias (bioucas@lx.it.pt), November, 2011
%
%%   

close all
clear all
verbose = 'off';

%% load data set  
load  Rterrain

[B,n] = size(Y);

% estimate noise
[w Rw]=estNoise(Y);

% p == subspace dimension
% Ek == subspace orthogonal basis
[p,Ek] = hysime(Y,w,Rw,verbose);

% project onto the signal subspace
Yp = Ek*Ek'*Y; 
clear Y;


% orthogonal projection
[Yorth,Uorth,my angles_orth, scales_orth] = affineProj(Yp,p,'proj_type','orth');
% dpft projection
[Ydpft,Udpft,my angles_dpft, scales_dpft] = affineProj(Yp,p,'proj_type','dpft'); %, ...
%                                                'u', Uorth(:,p)/(Uorth(:,p)'*my));

%  angle norm scattergram 
                                        
norm_yp = sqrt(sum(Yp.^2));
figure(1)
[angles_orth_descend index] = sort(angles_orth,'descend');
plot(norm_yp(index(1:10000)),angles_orth_descend(1:10000), '.' )
%title('Orthogonal proj: angle-norm scattergram','FontSize',16)
xlabel('||y||','FontSize',16)
ylabel('angle(y,yp) [degrees]','FontSize',16)
set(gca,'FontSize',16);
text(0.3,10, 'Rcuprite data set','FontSize',16)


norm_ydpft = sqrt(sum(Ydpft.^2));
figure(2)
plot(norm_yp(1:10:end), norm_ydpft(1:10:end),'.' )
title('Perspective proj: norm(y)-norm(yp) scattergram','FontSize',16)
xlabel('norm(y)','FontSize',16)
ylabel('norm(yp)','FontSize',16)
set(gca,'FontSize',16);
text(24,120, 'Rcuprite data set','FontSize',16)


