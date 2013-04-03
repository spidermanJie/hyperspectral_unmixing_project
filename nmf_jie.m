function [W,H]=nmf_jie(X,K,maxiter,speak,WW,TA,Data)
%
% NMF using alternating least squares with obtimal brain surgeon.
%
% INPUT:
% X (N,M) : N (dimensionallity) x M (samples) non negative input matrix
% K       : Number of components
% maxiter : Maximum number of iterations to run
% speak   : prints iteration count and changes in connectivity matrix 
%           elements unless speak is 0
%
% OUTPUT:
% W       : N x K matrix
% H       : K x M matrix
%
% Jie Feng
% Mathematics
% University of California,Irvine
% jfeng1115@gmail.com
% 2012/10/12

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User adjustable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print_iter = 50; % iterations between print on screen and convergence test
obs_steps  = 15; % number of OBS steps to run before truncation
TA = TA./(repmat(sum(TA),size(TA,1),1));
WW = WW./(repmat(sum(WW),size(WW,1),1));
Data = Data./(repmat(sum(Data),size(Data,1),1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for negative values in X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if min(min(X)) < 0
    error('Input matrix elements can not be negative');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize random W and H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kk = 3 ;
[D,N]=size(X);
W=rand(D,K);
H=rand(K,N);

% use W*H to test for convergence
Xr_old = W*H;
 dist = zeros (K,N);
 lamda = -0.005; 



for n=1:maxiter,
    W = ((pinv(H*H')*H)*X')';
    %%% find nearest kk %%%

    for ss = 1:K 
        dist(ss,:) = sum((repmat(W(:,ss),1,N) - X).^2); 
        [y,ind ] = min(dist(ss,:)); 
        dW(:,ss) = lamda * (W(:,ss)-y);
        
    end
    W = W+dW; 
    %%% OSB END %%%
    W = (W>0).*W; % truncate negative elements
    W = W./repmat(sum(W),D,1); % normalize columns to unit length

    %%%%%%%%%%%%%%%% 
    % H update
    %%%%%%%%%%%%%%%%
    H=(W*pinv(W'*W))'*X;
    %%% OSB %%%
    
    %%% END OBS %%%
    H=H.*(H>0); % truncate negative elements
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % print to screen
    %%%%%%%%%%%%%%%%%%%%%%%
    if (rem(n,print_iter)==0) & speak,
        W1 = W./(repmat(sum(W),size(W,1),1)); 
    
        figure(); 
        scatter3(W1(1,:),W1(2,:),W1(3,:),'b.'); 
         hold on; 
         scatter3(WW(1,:),WW(2,:),WW(3,:),'r.'); 
         hold on; 
         scatter3(TA(1,:),TA(2,:),TA(3,:),'g.'); 
          hold on; 
         scatter3(Data(1,:),Data(2,:),Data(3,:),'m.');
         hold on;
         legend('new algorithm','NMF','True','data vectors')
        Xr = W*H;
        diff = sum(sum(abs(Xr_old-Xr)));
        Xr_old = Xr;
        eucl_dist = nmf_euclidean_dist(X,W*H);
        errorx = mean(mean(abs(X-W*H)))/mean(mean(X));
        disp(['Iter = ',int2str(n),...
            ', relative error = ',num2str(errorx),...
            ', diff = ', num2str(diff),...
            ', eucl dist ' num2str(eucl_dist)])
        if errorx < 10^(-5), break, end
    end
end