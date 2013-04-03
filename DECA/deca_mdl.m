function [varargout]=deca_mdl(varargin)
%
%% -------------- DECA Sintax -------------- 
%
% Syntax: [M_deca,S_deca] = deca_mdl(y);
%         [M_deca,S_deca] = deca_mdl(y,Deca_Parameters);
%         [M_deca,S_deca,Up,parameters_estimates,log_likelihood_evolution] = deca_mdl(y,Deca_Parameters);
%
%% --------------  Required inputs  -------------- 
%
% y - Observed dataset [L x N] matrix,
%     where L is the number of channels (bands),
%     N is the number of pixels (Lines x Columns).
%     Each pixel is a linear mixture of p endmembers
%     signatures, i.e.,
%
%              y = M*s + noise,
%
%     where M is an [L x p] matrix with endmembers 
%     signatures in each column, 
%     s is an [p x N] matrix with the endmembers 
%     abundance fractions, which are assumed to follow 
%     a mixture of Dirichlet distributions.
%
%% -------------- Optional inputs -------------- 
%
% Deca_Parameters - Structure with optional parameters
%                   (the structure can be given with all or some parameters)
% Deca_Parameters = struct('Endmembers',3, ...
%                          'Number_mode_max',5, ...
%                          'Number_mode_min',1, ...
%                          'Tol_th',1e-5, ...
%                          'Number_iteration_max',400, ...
%                          'Initial_Endmembers_signatures',M_ini, ...
%                          'Initial_mode_weights', e_ini ...
%                          'Initial_Dirichlet_Parameters',a_ini, ...
%                          'Projection_matrix',Up, ...
%                          'Display_figure','on', ...
%                          'Image_Dimensions',[N_LINES N_COLUMNS N_BANDS],...
%                          'Verbose','on');
% where:
% 'Endmembers'           - Number of endmembers [positive integer]
% 'Number_mode_max'      - Maximum number of Dirichlet modes [positive integer]
% 'Number_mode_min'      - Minimum number of Dirichlet modes 
%                          (smaller than Kmod_max) [positive integer]
% 'Tol_th'               - Tolerance for the termination test (relative
%                          variation of the log-likelihood) [double]
% 'Number_iteration_max' - Maximum number of iterations [positive integer]
% 'Initial_Endmembers_signatures'
%                        - Initial estimate for the endmembers signatures
%                          [L x p] matrix 
%                          (optional parameter, if not given the initial estimate is given by SISAL)
% 'Initial_mode_weights'
%                        - Initial estimate for the mode weights 
%                          the sum must be equal to one
%                          [1 x Kmod_max] vector (optional parameter)
% 'Initial_Dirichlet_Paramenters'
%                        - Initial estimate for the Dirichlet parameters
%                          [p x Kmod_max] matrix (optional parameter)
% 'Projection_matrix'    - Dimensionality reduction [Lxp] matrix
%                          spanning the signal subspace
% 'Display_figure'       - Option to show figures while running DECA [string]
%                          'on' or 'off' default is 'on'
% 'Image_Dimensions'     - vector with the number of lines, columns, and
%                          bands of the dataset Y (projected or unprojected)
%                          [N_LINES N_COLUMNS N_BANDS]
% 'Verbose'              - Option to display information [string]
%                          'on' or 'off' default is 'on'
%
%% -------------- Outputs Parameters -------------- 
%
%  M_deca - Endmembers signatures estimates ([L x p] matrix)
%  S_deca - Abundance fractions estimates ([p x N] matrix)
%
%% -------------- Aditional Outputs -------------- 
% 
% Aditional output Up - Projection_matrix
%
% Output "parameters_estimates" show  
% the estimated number of modes, their weights, 
% their parameters, and the posterior probability.
%
% parameters_estimates = struct('Number_of_modes',k_opt,...
%                               'Dirichlet_parameters',a_opt,...
%                               'Dirichlet_weigths',e_opt,...
%                               'Post_probability',beta_opt);
%
% Output "log_likelihood_evolution" show the 
% evolution of the log_likelihood function, 
% the evolution of the Dirichlet weights,
% the evolution of the emdmembers signatures, and 
% the transitions where a mode is anihilated.
%
% log_likelihood_evolution = struct('log_likelihood',L_global,...
%                                   'weigths',e_global,...
%                                   'endmember_matrix',M_global,...
%                                   'transitions',transitions);
%  
%
%% -------------- Brief DECA Description -------------- 
%
% DECA is an usupervised  method to unmix highly linear 
% mixed hyperspectral datasets, in which the simplex of 
% minimum  volume,  usually  estimated  by  the  purely 
% geometric based  algorithms, is far way from the true 
% simplex associated with the endmembers.
% DECA uses a mixture of  Dirichlet  densities as prior 
% for the abundance  fractions, which  allows to  model 
% complex  distributions in which the mass  probability 
% is scattered by a number of clusters inside the simplex. 
% Furthermore,  the  Dirichlet  density   automatically 
% enforces  non-negativity and  sum-to-one  constraints
% on the abundance fractions.
% DECA consists in a cyclic minimization algorithm: 
% 1) the  number  of  Dirichlet  modes  are  adaptively
% inferred based on the minimum description length (MDL)
% principle;
% 2) a  generalized   expectation   maximization  (GEM)
% algorithm infers the model parameters;
% 3) a sequence of augmented Lagrangian based optimizations
% are used to compute the signatures of the endmembers.
%
% More details in:
%
% José M. P. Nascimento and José M. Bioucas-Dias,
% "Hyperspectral Unmixing based on Mixtures of Dirichlet Components"
% IEEE Transaction on Geoscience and Remote Sensing
% Vol. 50, Nº 3, pp. 863-878, 2012
%
%% -------------- Copyright -------------- 
%
% Copyright (2011):        
% José M.P. Nascimento (zen@isel.pt)
% & 
% José Bioucas-Dias (bioucas@lx.it.pt)
%
% Created: January 2011
% Latest Revision: -

%
% DECA is distributed under the terms of
% the GNU General Public License 2.0.
%
% Permission to use, copy, modify, and  distribute this 
% software for any purpose without fee is hereby granted,
% provided that  this entire  notice is included in all 
% copies of any software which is or includes a copy or 
% modification of this software and in all copies of the
% supporting documentation for such software.
% This  software is  being provided "as is", without any
% express or implied warranty. In particular, the authors
% do not make any representation or warranty of any kind 
% concerning the merchantability of this software or its 
% fitness for any particular purpose."
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --------------  Read input parameters  -------------- 
[Y,N_LINES,N_COLUMNS,L,N,p,U,Kmod_max,Kmod_min,A_est,e_hat,a_hat,Tol_th,Max_iter,vefig,verb] ...
    = read_deca_params(varargin);
%% -------------- looking for output parameters -------------- 
if nargout > 5, error('Too many output parameters'); end;
if (nargout < 2) & verb, fprintf(1,'Warning: Output only the endmembers signatures.\n');end;
if nargout == 5, rec_evolution=1; else rec_evolution=0;end

%% Estimate The signal subspace projection matrix (if it is not given)
if p~=L, % if p==L it means that the dataset is already projected onto the subspace
   if any(isnan(U(:)))
      if verb, fprintf(1,'Using Hysime to estimate the signal subspace.\n'); end;
      % Estimate noise matrix, assuming additive noise
      [Rw w]=HyperNoise(Y,'additive','off'); 
      % Estimate number of endmembers
      HySime_Parameters = struct('Noise',w,'Noise_Correlation_Matrix',Rw, ...
                                 'Display_figure','off','Verbose','off');
      [p_hat,U,out_struct]=HySime(Y,HySime_Parameters);   
      %E = out_struct.Projection_Matrix_Percent(:,1:p);

      % if the number of endmembers is not given, use the Hysime estimation
      if isnan(p)
         if verb,fprintf(1,'Using Hysime to estimate the number of endmembers.\n');end
         p = p_hat;
      else
         U = out_struct.Projection_Matrix_Percent(:,1:p);
      end
      clear w Rw HySime_Parameters out_struct p_hat  
   end
   %% -------------- Project the dataset to the signal subspace -------------- 
   Y = U'*Y;
end

%% -------------- Project the dataset to the affine set -------------- 
[Y,Up,my,sing_val] = dataProj(Y,p,'proj_type','affine');

%% -------------- Initial Endmember Signature Estimate by sisal -------------- 
if any(isnan(A_est(:)))
  if verb,fprintf(1,'Using SISAL to estimate initial endmembers signatures. \n');end
  [A_est, dummy, Z] = dusal(Y,p, 'dir_pars',ones(p,1),...
                                 'spherize', 'no',...
                                 'MM_ITERS',40,...
                                 'MU',1,...
                                 'verbose',0);
else
  if p~=L,  A_est = U'*A_est;  end
end
%% -------------- Show dataset and endmembers initial estimates -------------- 
if vefig
   figure(2)
     
     plot(Y(1,:),Y(2,:),'b.','markersize',3)
     hold on
     plot(A_est(1,1:p), A_est(2,1:p),'ro','markersize',8)
     hold off
end
%% -------------- Initial Dirichlet Parameters Estimate -------------- 
if isnan(e_hat)
   a_hat = 1+.5*rand(p,Kmod_max); 
   e_hat = ones(1,Kmod_max)/Kmod_max+1/(Kmod_max)*rand(1,Kmod_max);
   e_hat = e_hat/sum(e_hat);
end

%% -------------- DECA algorithm -------------- 
[A_opt,s_opt,parameters_estimates,log_likelihood_evolution] ...
    = deca_core(Y,N_LINES,N_COLUMNS,L,N,p,A_est,e_hat,a_hat,Kmod_max,Kmod_min,Tol_th,Max_iter,rec_evolution,vefig,verb);
if p~=L, 
    % back to the original subspace
    M_deca = U*A_opt; 
else
    M_deca = A_opt; 
end;
S_deca = s_opt;
%% -------------- output results -------------- 
varargout(1) = {M_deca};
varargout(2) = {S_deca};
varargout(3) = {U};
if nargout >= 4, varargout(4) = {parameters_estimates};end
if nargout == 5, varargout(5) = {log_likelihood_evolution};end
end
%% end of function [varargout]=deca_mdl(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INTERNAL FUNCTIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% deca_core - Core of DECA algorithm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A_opt,s_opt,parameters_estimates,log_likelihood_evolution] ...
        = deca_core(Y,N_LINES,N_COLUMNS,L,N,p,A_est,e_hat,a_hat,Kmod_max,Kmod_min,...
                    Tol_th,Max_iter,rec_evolution,vefig,verb)

        % default colors to draw figures
         ycolor =[0 .5 1; .5 0 1; 0 1 1;1 0 1; 0 1 .5; 0 .5 .5; .5 1 0; 1 0 .5;.5 .5 .5; 1 1 .5];
        % if rec_evolution==0 these variables do not change
        % no record of the parameters evolution
         L_global = [];
         e_global = [];
         transitions = [];
         M_global = [];
         if rec_evolution, 
            M_global=cat(3,M_global,A_est);
         end
         
         L_opt=+inf;
         Kmod =Kmod_max;
         % Run deca for each number of modes
         while Kmod >= Kmod_min
         
             L_old = -inf; 
             iter=0;
             stop_while = 0;
  
                    
             while ~stop_while  & (iter < Max_iter) 
                   iter = iter + 1;   
                   % estimate the abundance fractions
                   s_hat=pinv(A_est)*Y;
                   % preventing that all samples are inside the simplex
                   s_hat= s_hat.*(s_hat>=0) + 1e-6;
                   % estimate the Dirichlet parameters
                   [a_hat,e_hat,beta,L_tot,e_tot,kill_mode]=MixDiparEst(s_hat,Kmod,Tol_th,Max_iter,a_hat,e_hat,0*verb,0*vefig);
                   
                   if Kmod==1, dir_pars = a_hat;
                   else dir_pars = a_hat * beta';
                   end
                   % estimate the endmembers matrix 
                   [A_est] = dusal(Y,p, 'dir_pars', dir_pars ,'spherize', 'no','MM_ITERS',5,'MU',1 ,'M0',A_est,'verbose',0);
         
                   % Annihilate mode with zero weigth (e_hat)
                   % Here we could annihilate more than one mode
                   ind = 0;
                   if find(kill_mode),
                      ind = kill_mode(1) ;
                   else 
                      if find(e_hat < 1e-3);
                         [val ind]=min(e_hat);
                      end
                   end
                   if ind
                      e_hat=e_hat(find((1:Kmod)~=ind));
                      a_hat=a_hat(:,find((1:Kmod)~=ind));
                      if rec_evolution, transitions = [transitions; Kmod size(L_global,1)];end
                      Kmod = Kmod -1;
                   end
                   % set stop condition
                   L_new = L_tot(end) + N*log(abs(det(A_est)));
                   L_variat = (L_new - L_old)/L_old;
                   stop_while = (abs((L_old - L_new)/L_old) < Tol_th);
                   L_old = L_new;

                   if rec_evolution, 
                      L_global = [L_global; L_new];
                      e_global = [e_global; [sort(e_hat,2,'descend')  nan*ones(size(e_hat(end,:),1),Kmod_max-Kmod)]];           
                      M_global = cat(3,M_global,A_est); 
                   end
                   % show the evolution 
                   if vefig
                      th = .99;
                      figure(2*Kmod_max+Kmod)
                         idx=find(prod(double(beta<=th)'));
                          hold on
                          plot(Y(1,idx),Y(2,idx),'k.','markersize',3)
                          for k = 1:Kmod,
                              idx = find(beta(:,k)>th);
                              if ~isempty(idx)
                                 plot(Y(1,idx),Y(2,idx),'.','Color',ycolor(k,:),'markersize',3)
                              end
                          end
                          plot(A_est(1,1:p), A_est(2,1:p),'.', 'Color',[0.8 0.5 0.3])
                          hold off
                          drawnow
                          figure(Kmod_max+Kmod)
                        subplot(211);plot(L_global);
                        title('likelihood');
                        subplot(212);plot(e_global);
                        title('the mixing weights' ); 
                      
       
                   end
                   if (verb & mod(iter,20)==0)
                       iter
                       Kmod
                       [e_hat;nan*ones(size(e_hat)); a_hat] 
                   end
             
                   % store the best parameters
                   if L_new < L_opt
                      k_opt = Kmod;
                      L_opt = L_new;
                      a_opt = a_hat;
                      e_opt = e_hat;
                      beta_opt = beta;
                
                      A_opt = A_est;
                      s_opt = s_hat;
                   end
             end 
             %Annihilate mode with smallest weigth e_hat
             [val ind]=min(e_hat);
             e_hat=e_hat(find((1:Kmod)~=ind));
             a_hat=a_hat(:,find((1:Kmod)~=ind));
             if rec_evolution, 
                 transitions = [transitions; Kmod size(L_global,1)];
             end
             Kmod = Kmod -1;
         
         end
         % show the final endmembers matrix estimate 
         if vefig
            figure(2*Kmod_max+Kmod)
              
              hold on;
              plot(A_est(1,[1:p 1]),A_est(2,[1:p 1]),'k*--','Color','r')
              hold off
         end
         % output variables
         parameters_estimates = struct('Number_of_modes',k_opt,...
                                 'Dirichlet_parameters',a_opt,...
                                 'Dirichlet_weigths',e_opt,...
                                 'Post_probability',beta_opt);
         log_likelihood_evolution = struct('log_likelihood',L_global,...
                                           'weigths',e_global,...
                                           'endmember_matrix',M_global,...
                                           'transitions',transitions);
end % function deca_core

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% read_deca_params - read all input parameters 
%                    and set default values if needed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Y,N_LINES,N_COLUMNS,L,N,p,U,Kmod_max,Kmod_min,M_ini,e_ini,a_ini,Tol_th,Max_iter,vefig,verb]=read_deca_params(varargin)

         % default parameters 
            p = nan; % estimation is needed
            U = nan; % estimation is needed
            M_ini = nan; % estimation is needed
            e_ini = nan; % estimation is needed
            a_ini = nan; % estimation is needed
            Kmod_e = nan; % auxiliary variable
            Kmod_max = 5;
            Kmod_min = 1;
            Tol_th = 1e-6;
            Max_iter = 500;
            vefig = 1; 
            verb = 1;

         input_data = varargin{1}; % varargin is a cell with 2 cells inside
         n_input = numel(input_data);
         if n_input < 1,
            error('Input observed data matrix is needed');
         end
         % reading required input matrix (observed data)
         Y = input_data{1};

         if isempty(Y) | ~isa(Y,'double') | any(isnan(Y(:))) | ndims(Y)~=2, 
            error('Wrong data matrix: Y should be an L (channels) x N (pixels) matrix, with double values.\n');
         end
         [L N] = size(Y);  % L number of bands (channels)
                           % N number of pixels (Lines x Columns)
         N_LINES = N; N_COLUMNS = 1; % default values (if not given)
         if L > N, 
            fprintf(1,'Warning: Aparently the number of observations is smaller than the number of channels. \n')
         end
         if n_input > 2,
            fprintf(1,'Warning: Too many input parameters. Using default parameters instead. \n')
         end
         if (n_input == 1)
            fprintf(1,'Using default parameters. \n');
         end
         % checking for optional parameters struct
         if (n_input == 2)
            Deca_Parameters = input_data{2};
            if ~isstruct(Deca_Parameters), 
               fprintf(1,'Warning: Invalid parameters struct. Using default parameters instead. \n');
            else
                % 'try' do not need 'catch', default parameters are already defined
                try % reading verbose first 
                    Verbose_in = lower(Deca_Parameters.Verbose); 
                    if ~any(strcmp(Verbose_in,{'on','off'}))
                       fprintf(1,'Warning: wrong Verbose option. Using default option. \n');
                    else verb = strcmp(Verbose_in,'on');
                    end
                end
                try 
                   Dim=Deca_Parameters.Image_Dimensions;
                   if (isnan(Dim(:))| numel(Dim)~=3 | sum(Dim<=0) | sum(rem(Dim,1)~=0))
                      if verb,fprintf(1,'Warning: Wrong Image_Dimension parameter. \n');end
                   else
                      N_LINES = Dim(1);
                      N_COLUMNS = Dim(2); 
                      if (N ~= N_LINES * N_COLUMNS ) % | L ~= Dim(3) % not included 
                         if verb,fprintf(1,'Warning: Image_Dimension parameter is not compatible with data matrix dimensions . \n');end                          
                         N_LINES = N; N_COLUMNS = 1; 
                      end
                   end
                end
                try 
                   p_in=Deca_Parameters.Endmembers;
                   if (p_in<0 | p_in>L | rem(p_in,1)~=0),  
                      if verb,fprintf(1,'Warning: ENDMEMBER parameter must be positive integer smaller than the number of channels. \n');end
                   else
                      p=p_in;
                   end
                   % if p is not given as input parameter, 
                   % it can be extracted from the Projection_matrix dimensions
                   % or from the Initial_Endmembers_signatures matrix
                   % or from the initial Dirichlet parameters
                end
                try % Assume that p is already verified!
                    U_in = Deca_Parameters.Projection_matrix;
                    [Lu pu]=size(U_in);
                    if any(isnan(U_in(:))) | ~isa(U_in,'double'),   
                       if verb,fprintf(1,'Wrong Projection matrix. \n');end
                    else
                       if L~=Lu
                          if verb,fprintf(1,'Projection matrix is not consistent with data matrix. \n');end
                       end
                       if isnan(p),
                          p = pu; % defining p by Projection matrix dimensions
                          if verb,fprintf(1,'Assuming the number of endmembers (p = %d) by Projection matrix inspection. \n',p);end
                          U=U_in;
                       elseif p~=pu,
                          if verb,fprintf(1,'Projection matrix is not consistent with number of endmembers. \n');end
                       else
                          U=U_in;
                       end
                    end
                end
                try % Assume that p is already verified!
                    M_in = Deca_Parameters.Initial_Endmembers_signatures;
                    [Lm pm]=size(M_in);
                    if any(isnan(M_in(:))) | ~isa(M_in,'double') | L~=Lm,   
                       if verb,fprintf(1,'Wrong initial endmembers estimate. \n');end
                    else
                       if isnan(p),
                          p = pm; % defining p by initial endmembers estimate
                          if verb,fprintf(1,'Assuming the number of endmembers (p = %d) by initial endmembers signatures. \n',p);end
                          M_ini=M_in;
                       elseif p~=pm,
                          if verb,fprintf(1,'Initial endmembers signatures not consistent with number of endmembers. \n');end
                       else
                          M_ini=M_in;
                       end
                    end
                end
                try 
                    e_in = Deca_Parameters.Initial_mode_weights;
                    [v ke]=size(e_in); 
                    if any(isnan(e_in(:))) | ~isa(e_in,'double') | v~=1 | abs(sum(e_in)-1)>1e-6
                       if verb,fprintf(1,'Wrong initial estimate mode weights. Using default estimates.\n');end
                    else
                       e_ini = e_in;
                       Kmod_e = ke; % number max of modes is assumed (verify later)
                    end
                end
                try % Assume that p and e_ini are already verified!
                    a_in = Deca_Parameters.Initial_Dirichlet_Parameters;
                    [pa ka]=size(a_in);
                    if any(isnan(a_in(:))) | ~isa(a_in,'double') | Kmod_e~=ka 
                       if verb,fprintf(1,'Wrong initial Dirichlet parameters. Using default estimates.\n');end
                    else
                       if isnan(p),
                          p = pa; % defining p by the parameters dimensions
                          if verb,fprintf(1,'Assuming the number of endmembers (p = %d) by initial Dirichlet parameters. \n',p);end
                          a_ini = a_in;
                       elseif p~=pa,
                          if verb,fprintf(1,'Initial Dirichlet parameters are not consistent with number of endmembers. \n');end
                       else
                          a_ini = a_in;
                       end
                    end
                end
                % additional test
                if any(isnan(e_ini(:))) | any(isnan(a_ini(:)))
                   kmod_e = nan;e_ini=nan;a_ini=nan;
                end
                %
                try % assuming that e_ini and a_ini are already verified
                    Kmod_max_in = Deca_Parameters.Number_mode_max; 
                    if rem(Kmod_max_in,1)~=0 | Kmod_max_in<=0, 
                       if verb,fprintf(1,'Warning: maximum number of modes must be a positive integer. Using default value: %d \n',Kmod_max);end;
                    elseif isnan(kmod_e) | Kmod_max_in==Kmod_e
                           Kmod_max = Kmod_max_in;
                    else
                           Kmod_max = Kmod_e;
                           if verb,fprintf(1,'Warning: maximum number of modes not consistent with the initial parameters. Assuming value: %d \n',Kmod_max);end;
                    end
                end
                try 
                    Kmod_min_in = Deca_Parameters.Number_mode_min; 
                    if rem(Kmod_min_in,1)~=0 | Kmod_min_in<=0, 
                       if verb,fprintf(1,'Warning: minimum number of modes must be a positive integer. Using default value: %d \n',Kmod_min);end
                    else Kmod_min = Kmod_min_in;
                    end
                end
                % additional test
                if Kmod_max < Kmod_min, 
                   if verb,fprintf(1,'Warning: Apparently the mode range is invalid: Number_mode_min < Number_mode_max\n');end
                   dummy = Kmod_max; Kmod_max = Kmod_min; Kmod_min = dummy;
                end
                %
                try 
                    Tol_th_in = Deca_Parameters.Tol_th; 
                    if ~isa(Tol_th_in,'double') | Tol_th_in<=0, 
                       if verb,fprintf(1,'Warning: wrong threshold value, Using default value: %d \n', Tol_th );end
                    else Tol_th = Tol_th_in;
                    end
                end
                try 
                    Max_iter_in = Deca_Parameters.Number_iteration_max; 
                    if rem(Max_iter_in,1)~=0 | Max_iter_in<=0, 
                       if verb,fprintf(1,'Warning: number maximum of iterations must be a positive integer. Using default value: %d \n',Max_iter);end
                    else Max_iter = Max_iter_in;
                    end
                end
                try 
                    Ver_fig_in = lower(Deca_Parameters.Display_figure); 
                    if ~any(strcmp(Ver_fig_in,{'on','off'}))
                       if verb,fprintf(1,'Warning: wrong Display_figure option. Using default option. \n');end
                    else vefig = strcmp(Ver_fig_in,'on');
                    end
                end
            end
         end
end % function read_deca_params

