%% Name: sparse_regression
%
%  Generate the parse reconstruction results plotted in Figure 19.
%  
%  Data:  simulated data set generated from the USGS library. 
%  
%  Top: Signal to reconstruction error (SRE) as a function of the
%  number of active materials. 
%  
%  Bottom: Number of incorrect selected material as a function of 
%  the number of active materials.
%
%
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Jose M. Bioucas-Dias, (bioucas@lx.it.pt), February, 2012)




%% begining of the simulation

clear all,
close all

% initializa random number generator
rand('seed',31416);
randn('seed',31416 );



% USGS_pruned_3_deg.mat is a pruned version  of the USGS library where
% any pair of signatures has angle not smaller than 3 degrees
%
load USGS_pruned_3_deg.mat
A = B;
[L,n] = size(A);


% reorder the signatues  by decrasing angles with any other signatures,
% that reorder the library by 
%
%   min_angle(i) =    min     angle(A(i,:),A(j,:))
%                  j ~= i
%

% find angles
angles = acos(((A'*A)./sqrt(sum(A.^2)'*sum(A.^2)))-1e-10);
angles = angles + 10*eye(n);
min_angle = min(angles);
[min_angle index_ang] = sort(min_angle, 'descend');
% library ordered by decreasing angles
AO = A(:,index_ang);

%
% AO(:,j) has an angle not smaller then min_angle(i) with any other
% signature

%%  set the parameters
N = 1000;   % number of pixels


for demo=1:6
    if demo == 1
        % demo 1
        SNR = 1000;     % no noise
        p0 = 1;         % select materials with  min_angle  >= 7 deg
        delta = 4;      % increment of the number of materials 
        lambda = 0;     % CLS
        tol = 1e-8;
        
    elseif demo == 2
        % demo 2
        SNR = 1000;    % no noise
        p0 = 200;      % select materials with  min_angle  <=  4 deg
        delta = 4;
        lambda = 0;    %CLS
        tol = 1e-8;
        
        
        % demo 3
    elseif demo == 3
        SNR = 25;
        p0 = 1;
        delta = 1;
        lambda = 5e-5;  % CSR
        tol = 1e-4;
        
        
    elseif demo == 4
        % demo 4
        SNR = 25;
        p0 = 1;
        delta = 1;
        lambda = 0;     % CLS
        tol = 1e-4;
        
    elseif demo == 5
        % demo 5
        SNR = 25;
        p0 = 200;
        delta = 1;
        lambda = 5e-2;  % CSR
        tol = 1e-4;
        
    elseif demo == 6
        % demo 6
        SNR = 25;
        p0 = 200;
        delta = 1;
        lambda = 0;     %CLS
        tol = 1e-4;
    end
    
    
   
    
    SHAPE_PARAMETER = 2;    % higly mixed data set
    MAX_PURIRY = 1;         % do not theshold abundances
    OUTLIERS   = 0;         % no outliers
    PURE_PIXELS = 'no';     % no pixels
    
    
    
    
    for i=1:10
        p = delta*i;
        supp = p0:p0+p-1;   % select signatures
        M = AO(:,supp);
        
        % generate the data
        [Y,x,noise] = spectMixGen(M,N,'Source_pdf', 'Diri_id','pdf_pars',SHAPE_PARAMETER,...
            'max_purity',MAX_PURIRY*ones(1,p),'no_outliers',OUTLIERS, ...
            'pure_pixels', PURE_PIXELS,'violation_extremes',[1,1.2],'snr', SNR, ...
            'noise_shape','uniform'); % % % % % % % % % % % % % % % % % % % % % % % % %
        
        
        % do sparse regression
        [xhat] = sunsal(AO,Y, 'POSITIVITY', 'yes', 'verbose', 'yes', ...
            'ADDONE','no', 'lambda',lambda,'TOL', tol);
        
        % detect the support
        [aux index] = sort(sum(xhat,2),'descend');
        
        aux  = length(setdiff(index(1:p), supp))
        err_supp(demo,i) = aux;
        
        % compute the SRE 
        aux = 20*log10(norm(x,'fro')/norm(xhat(supp,:)-x,'fro'))
        SRE(demo,i) = aux;
        
        
        
    end
    
end

%save SRE_err_supp SRE err_supp

delta = 4 
figure(1)
plot(delta:delta:10*delta, SRE(1,:), delta:delta:10*delta, SRE(2,:), ...
       'LineWidth', 2)
   
xlabel('||x||_0 (number of active columns of A)','FontSize',16)   
ylabel('SRE (dB)','FontSize',16)   
legend('\theta_{min}(A) \geq 7^\circ', '\theta_{min}(A) \leq 4^\circ')
title('Signal to Reconstruction Error', 'FontSize',16)
text(5,35, 'SNR = \infty','FontSize',16)
text(5,19, '\lambda = 0','FontSize',16)
set(gca,'FontSize',16)


   
delta = 1
figure(2)
plot(delta:delta:10*delta, SRE(3,:),'b-', delta:delta:10*delta, SRE(4,:), 'b--', ...
      delta:delta:10*delta, SRE(5,:),'r-', delta:delta:10*delta, SRE(6,:),'r--','LineWidth', 2)
   
xlabel('||x||_0 (number of active materials)','FontSize',16)   
ylabel('SRE (dB)','FontSize',16)   
legend('\theta_{min}(A) \geq 7^\circ, \lambda = 5\times 10^{-4}', '\theta_{min}(A) \geq 7^\circ, \lambda = 0', ...
       '\theta_{min}(A) \leq 4^\circ, \lambda = 5\times 10^{-2}', '\theta_{min}(A) \leq 4^\circ, \lambda = 0')
title('Signal to Reconstruction Error', 'FontSize',16)
text(0.5,20, 'SNR = 25 dB','FontSize',16)
set(gca,'FontSize',16)

  

delta = 1
figure(3)
plot(delta:delta:10*delta, err_supp(3,:), 'b-', delta:delta:10*delta, err_supp(4,:), 'b--',...
      delta:delta:10*delta, err_supp(5,:),'r-', delta:delta:10*delta, err_supp(6,:),'r--', 'LineWidth', 2)
   

xlabel('||x||_0 (number of active materials)','FontSize',16)   
title('Number of incorrect materials','FontSize',16)   
legend('\theta_{min}(A) \geq 7^\circ, \lambda = 5\times 10^{-4}', '\theta_{min}(A) \geq 7^\circ, \lambda = 0', ...
       '\theta_{min}(A) \leq 4^\circ, \lambda = 5\times 10^{-2}', '\theta_{min}(A) \leq 4^\circ, \lambda = 0')
text(0.4,45, 'SNR = 25 dB','FontSize',16)
set(gca,'FontSize',16)

  
  
