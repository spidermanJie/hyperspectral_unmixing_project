clc;
clear all;
close all; 
addpath /Users/jfeng1115/Documents/MATLAB/Demix
m = 3; 
p = 3; 
N =100; 
Talpha = [39,7,3;3,8,40;7,35,4]' ; % the True parameter of Dir distribution
Teps = randi([2,20],1,m);% the True parameter for the mixture weights 
Teps = Teps/sum(Teps); 
r = mnrnd(1,Teps,N) ;
Source = zeros(N,p);
sum2 = sum(r); 
for i =1:m
    ind = find(r(:,i)); 
    S1 = sample_dirichlet(Talpha(i,:),sum2(i)); 
    Source(ind,:) = S1; 
    
end
Source =formatData(Source);
alpha = 4:1:6;
NN = length(alpha);
T = 100;
ss = rand(NN,T);
for n = 1 : NN
    ss(n,:) = gamrnd(alpha(n),1,1,T);    % these are Gamma sources
end
S = Source'*mean(mean(ss))/mean(mean(Source));
S = [S ss];

K = 3; 

%S = [S,[2.9,0.14,0.22]',[0.14,3.8,0.22]',[0.29,0.14,4]'];
% the matrix W
%f is the number of data signal 

%generate ture TW and TA
f = 3; 
p = 3;
TA = abs(zeros(f,p)); 
while (det(TA)<0.1)
TA = abs(rand(f,p));
end
%TW = inv(TA); 
figure(1);

subplot(3,1,1), plot(S(1,:)); title ('source signal1');
subplot(3,1,2), plot(S(2,:)); title ('source signal2');
subplot(3,1,3), plot(S(3,:)); title ('source signal3');
title('source')

%subplot(2,2,4), plot(S(4,:)); 
%generate Data
Data = TA*S;
TA
S
[W,H] = nmf(Data,K,'als',200);
[W1,H1]=nmf_jie(Data,K,200,1,W,TA,Data);
W = W./(repmat(sum(W),size(W,1),1)); 
W1 = W1./(repmat(sum(W1),size(W1,1),1)); 
TA = TA./(repmat(sum(TA),size(TA,1),1));
Data = Data./(repmat(sum(Data),size(Data,1),1));
figure();
title('NMF');
subplot(3,1,1), plot(H(3,:)); 
subplot(3,1,2), plot(H(2,:)); 
subplot(3,1,3), plot(H(1,:)); 

figure();
title('new algorithm');
subplot(3,1,1), plot(H1(3,:)); title('recovered signal1');
subplot(3,1,2), plot(H1(1,:)); title('recovered signal2');
subplot(3,1,3), plot(H1(2,:)); title('recovered signal3');
%subplot(2,2,4), plot(H(4,:)); 

figure(); 
scatter3(W1(1,:),W1(2,:),W1(3,:),'b.'); 
hold on; 
scatter3(W(1,:),W(2,:),W(3,:),'r.'); 
hold on; 
scatter3(TA(1,:),TA(2,:),TA(3,:),'g.'); 
hold on; 
scatter3(Data(1,:),Data(2,:),Data(3,:),'m.');
legend('new algorithm','NMF','True','data vectors')
TA
W1
W

hold off; 




