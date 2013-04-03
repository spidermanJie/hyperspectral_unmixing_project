function [a_hat,e_hat,b,L_tot,e_hat_tot,kill_mode]=MixDiparEst(s_hat,Q,th,niter_max,a_hat,e_hat,verbose,ver_fig)

%% -------------- MixDiparEst sintax -------------- 
%
% Syntax:
%        [a_hat,e_hat,b,L_tot,e_hat_tot]
%        = MixDiparEst(s,Q,th,niter_max,a_ini,e_ini,ver_fig,verbose)
%
% Mixture of Dirichlet parameters estimation
% This code is part of the DECA method.
%
%
%% -------------- Input Parameters -------------- 
% 
% s - [p x N] matrix with abundance fractions (s.t. s_ij>0 and sum_i(s_ij)=1)
%
% Q - number of Dirichlet modes [double]
%
% th - Tolerance for the termination test (relative
%      variation of the log-likelihood) [double]
%
% niter_max - Maximum number of iterations [positive integer]
%
% a_ini - Initial estimate for the Dirichlet parameters
%         [p x Kmod_max] matrix (optional parameter)
%
% e_ini - Initial estimate for the mode weights 
%         the sum must be equal to one
%         [1 x Kmod_max] vector 
% 
% verbose - Option to display information [string]
%
%% %%%%%%%%%%%%%%%% Output Parameters %%%%%%%%%%%%%%%%
%
% a_hat     - estimated mode parameters 
% e_hat     - estimated mode weights
% b         - posterior probability
% L_tot     - evolution of the log_likelihood function
% e_hat_tot - evolution of the Dirichlet weights
%
%% ------------------------------------------ 
% More details in:
%
% José M. P. Nascimento and José M. Bioucas-Dias,
% "Hyperspectral Unmixing based on Mixtures of Dirichlet Components"
% IEEE Transaction on Geoscience and Remote Sensing
% Vol. 50, Nº 3, pp. 863-878, 2012
%
%  For any comments contact the authors
%% -------------- Copyright -------------- 
%
% Copyright (January, 2009):
%
%  José Nascimento (zen@isel.pt)
%  & 
%  José Bioucas-Dias (bioucas@lx.it.pt)
%
% Latest Revision: May 2011
%
% MixDiparEst is distributed under the terms of
% the GNU General Public License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ----------------------------------------------------------------------

[p N] = size(s_hat);

if verbose,
   n_out = sum((sum(s_hat)-1)>1e-5);
   fprintf(1,'number of observations out of the simplex: %d\n', n_out);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EM algorithm to estimate the
% parameters of Dirichlet mixtures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% allocate memory to store the loglikelihhod 
% and the mixing weights evolution
e_hat_tot = nan*zeros(niter_max,Q);
L_tot = nan*zeros(niter_max,1);

L_old = -inf;
iter=0;
stop_while_Lfunc = 0; 
stop_while_modezero = 0;
while ~stop_while_Lfunc  & (iter < niter_max) & ~stop_while_modezero
    iter = iter + 1;   
    % E-step
    % Equation (4)
    c = exp(gammaln(sum(a_hat+realmin))-sum(gammaln(a_hat+realmin)));
    kill_mode = find(c > 1e200);
    stop_while_modezero=any(kill_mode); %
    c_r = repmat(c,[N 1]);
    a_hat_r = reshape(repmat(a_hat,[N 1]),[p N Q]);
    s_hat_r = reshape(repmat(s_hat,[1 Q]),[p N Q]);
    if Q==1
        D = c_r .* (prod(s_hat_r.^(a_hat_r-1),1))';        
    else
        D = c_r .* squeeze(prod(s_hat_r.^(a_hat_r-1),1));
    end
    e_r = repmat(e_hat,[N 1]);
    eD = e_r.*D ;
    % Equation (8)
    b = eD./repmat(sum(eD,2),[1 Q]);
    % M-step
    % Equation (10)
    e_hat = mean(b);
    % Equation (11)
    b_r = reshape(repmat(reshape(b,[N*Q 1])',[p 1]),[p N Q]);
    mbx = squeeze(sum(b_r.*log(s_hat_r),2))./repmat(sum(b)+realmin,[p 1]);
    a_hat = psi_inv( repmat(psi(sum(a_hat)+realmin),[p 1]) + mbx  );
    % Equation (17)
    L_prob = - sum(log(sum(eD,2))); % the term -N*log(abs(det(W))) is constant
    L_pen = Q*(p+1)/2*(1+log(N/12)) + p/2*sum(log(e_hat+realmin));
    %%% Q*(p+1)/2*(log(N)) + 1/2*sum(log(e_hat+realmin)); very old version 
    L_new = L_prob + L_pen;
    % set stop condition
    stop_while_Lfunc = (abs((L_old-L_new)/L_old) < th);
    L_old = L_new;
    L_tot(iter,:) = L_new; 
    e_hat_tot(iter,:) = e_hat;
    % display result
    if verbose, [e_hat; nan*ones(1,Q); a_hat ], end;
    % show loglikelihood evolution 
    % and mixing mode weights 
    if ver_fig
       figure(3)
         subplot(211);plot(L_tot);
         subplot(212);plot(e_hat_tot)
       drawnow
    end
end % while
L_tot = L_tot(1:iter,:);
e_hat_tot = e_hat_tot(1:iter,:);
end % end of function []=MixDiparEst()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function param_est = est_dirichlet_param(x,verbose)

% Estimates the one-mode Dirichlet parameters using the moments method
% Receives a matrix x (n_vectors x n_observations)
% the sum of the vectors is 1 and they are Dirichlet distributed

if verbose,
   n_out = sum((sum(x)-1)>1e-5);
   fprintf(1,'number of observations out of the simplex: %d\n', n_out);
end

m = mean(x,2)';
v = var(x',1);

[vm i] = max(v);
alfa_0= m(i)*( 1-m(i) )/vm -1;	

param_est = m * alfa_0;

end
