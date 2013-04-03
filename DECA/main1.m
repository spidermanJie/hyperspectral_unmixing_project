%% -------------- MainDemoDECA -------------- 
%
% This demo illustrates how to run DECA   
% (DEpendent Component Analysis)
%
% See details in:
%
% José M. P. Nascimento and José M. Bioucas-Dias,
% "Hyperspectral Unmixing based on Mixtures of Dirichlet Components"
% IEEE Transaction on Geoscience and Remote Sensing
% Vol. 50, Nº 3, pp. 863-878, 2012
%
%% -------------- Brief DECA Description -------------- 
%3
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

clear all,
clc
close all
%% -------------- Simulation Options -------------- 
% option = 1 runs DECA on simulated data
% option = 2 runs DECA on Cuprite dataset

option = 1;

%% -------------- Simulation Parameters -------------- 

switch option

        
  case {1}   
        p = 3;  
        % number of endmembers
        mm = 3; %number of mixtures 
        N_LINES=50;                 % number of lines of the image
        N_COLUMNS=60;               % number of columns of the image 
        SNR = 2000;                  % signal-to-noise ratio (E ||y||^2/E ||n||^2) in dBs
        %SHAPE_PARAMETER = [4/10 6 25 9;
        %                   5/10 7 8 23
        %                   1/10 20 20  4];
        SHAPE_PARAMETER = [1/2 10 24 4;
                           1/5 8 14 34;
                           1-1/5-1/2 13 23 16];
                           
        Kmod_max = 3; 
        Kmod_min = 3; 
                       
       % Abundance fractions are Dirichlet distributed
       % SHAPE_PARAMETER determines their distribution over the simplex.
       %
       % Each line of SHAPE_PARAMETER correspond to each mode
       % first column contains the weight (SHAPE_PARAMETER(:,1)) 
       % the remaining columns (SHAPE_PARAMETER(i,2:p)) contain 
       % the parameters of a Dirichlet mode. 
       % where:
       %       0 < SHAPE_PARAMETER(i,1)<= 1 , 
       %       sum(SHAPE_PARAMETER(:,1)) = 1  , and
       %       SHAPE_PARAMETER(i,2:p) > 0, 
       % when SHAPE_PARAMETER are
       %   = 1   ->  uniform over the simplex
       %   > 1   ->  samples moves towards the center of the simplex, 
       %             corresponding to highly mixed materials.
       %   ]0,1[ ->  samples moves towards the facets of the simplex.

       %% -------------- Select p Signatures from USGS -------------- 
        load USGS_pruned_10_deg
        [L n_materiais]=size(M);
       
        sel_mat = randperm(n_materiais);
        set_2 = randperm(L);
        BANDS = set_2(1:mm);
        sel_mat = sel_mat(1:p);
        M = M(BANDS,sel_mat);
        wavlen = wavlen(BANDS,:);
        clear datalib names aux st n_materiais sel_mat 

       %% -------------- Generate the Hyperspectral Observations -------------- 
       % Observations are linear mixtures of the endmembers signatures plus noise

        N = N_LINES * N_COLUMNS;      % number of pixels
        [Y,x,noise,sigma,outliers,mode_length] = spectMixGen(M,N,'Source_pdf','Diri_mix','pdf_pars',SHAPE_PARAMETER,'snr',SNR);

        aux = cumsum([1; mode_length]);
        real_Kmod = length(mode_length);
        for i=1:real_Kmod
            eval(['Y' int2str(i) '=Y(:,aux(' int2str(i) '):aux(' int2str(i+1) ')-1);']);
        end
end % end of switch

%% -------------- Project the dataset to the signal subspace -------------- 
        Y = Y;
        m_true =M;
       

%% -------------- Project the dataset to the affine set -------------- 
         [Y,U_aff,my,sing_val] = dataProj(Y,p,'proj_type','affine');


%% -------------- Scatterplot -------------- 
      figure(2)
      if option==1
         bcolor =[0 .5 1; .5 0 1; 0 1 1;1 0 1; 0 1 .5; 0 .5 .5; .5 1 0; 1 0 .5;.5 .5 .5; 1 1 .5];
         %bcolor = [zeros(Kmod,1) linspace(0,1,Kmod)' ones(Kmod,1)];
         leg_str = [];
         for i=1:real_Kmod
             eval(['plot(Y' int2str(i) '(1,:),' 'Y' int2str(i) '(2,:),''.'',''Color'',bcolor(' int2str(i) ',:),''markersize'',5)']);
             hold on;
             leg_str=[leg_str '''data region' int2str(i) ''','];
         end
      end
   
      plot(m_true(1,[1:p 1]),m_true(2,[1:p 1]),'k*--','markersize',8,'Linewidth',2)
      hold on;
      title('True regions')
      
      for Kmod = Kmod_min:Kmod_max
        figure(Kmod_max*2+Kmod)
      if option==1
         bcolor =[0 .5 1; .5 0 1; 0 1 1;1 0 1; 0 1 .5; 0 .5 .5; .5 1 0; 1 0 .5;.5 .5 .5; 1 1 .5];
         %bcolor = [zeros(Kmod,1) linspace(0,1,Kmod)' ones(Kmod,1)];
         leg_str = [];
         for i=1:real_Kmod
             eval(['plot(Y' int2str(i) '(1,:),' 'Y' int2str(i) '(2,:),''.'',''Color'',bcolor(' int2str(i) ',:),''markersize'',5)']);
             hold on;
             leg_str=[leg_str '''data region' int2str(i) ''','];
         end
      end
      plot(m_true(1,[1:p 1]),m_true(2,[1:p 1]),'k*--','markersize',8,'Linewidth',2)
      hold on;
      title('True regions')
      end
      %eval( ['legend(' leg_str  '''Endmembers' ''')'])
%% -------------- plot endmembers signatures -------------- 
      %M_aux = nan*zeros(length(wavlen),p);
      %M_aux(BANDS,:) = M;
      M_aux=M;
      
      p_color=['b','m','r','g'];
      
      figure(1);
      subplot(212);
      hold on;
      for k=1:p
          plot(wavlen,M_aux(:,k),'-','Color',p_color(k));
      end
      axis([min(wavlen) max(wavlen) 0 1])
      hold off
      title('Signatures')
      xlabel('Wavelength (\mu m)')
      ylabel('Reflectance')
%% -------------- plot image band -------------- 
      
    


%% -------------- DECA parameters settings -------------- 

   Deca_Parameters = struct('Endmembers',p, ...                
                            'Number_mode_max',Kmod_max, ...
                            'Number_mode_min',Kmod_min, ...
                            'Tol_th',1e-7, ...
                            'Number_iteration_max',100, ...
                            ...'Projection_matrix',U, ...
                            ...'Initial_mode_weights', e_ini ...
                            ...'Initial_Dirichlet_Parameters',a_ini, ...
                            ...'Initial_Endmembers_signatures',M_ini, ...
                            'Display_figure','on', ...
                            'Image_Dimensions',[N_LINES N_COLUMNS p],...
                            'Verbose','on'...
                            );

%% -------------- DECA algorithm -------------- 
tic;                     
[M_deca,S_deca,Up,parameters_estimates,log_likelihood_evolution]=deca_mdl(Y,Deca_Parameters);
Tdeca=toc;
%% -------------- Show DECA results -------------- 

   fprintf(1,'DECA time: %f\n',Tdeca);
   fprintf(1,'Number of modes estimates: %d\n', parameters_estimates.Number_of_modes);
   for i=1:parameters_estimates.Number_of_modes
       fprintf(1,'mode %d weight: %.3f; ', i, parameters_estimates.Dirichlet_weigths(i));
       fprintf(1,'Dirichlet parameters estimates:[');
       fprintf(1,'%.3f ', parameters_estimates.Dirichlet_parameters(:,i));
       fprintf(1,']\n')
   end


   [Ld pd]=size(M_deca);
   if Ld==pd, % dataset projected before DECA
      m_deca = M_deca;
      M_deca = M_deca;
    else
      Y = Up'*Y;
      m_deca = Up'*M_deca;
   end

%% -------------- plot loglikelihhod evolution --------------     
   m_evolut = log_likelihood_evolution.endmember_matrix;
   k_opt = parameters_estimates.Number_of_modes;

   figure(1)
      subplot(211)
      plot(log_likelihood_evolution.log_likelihood) 
      title('Loglikelihood Evolution')
      xlabel('Number of iterations')
      ax=axis;
      for k=1:size(log_likelihood_evolution.transitions,1)
          line([log_likelihood_evolution.transitions(k,2) ...
                log_likelihood_evolution.transitions(k,2)],[ax(3) ax(4)],'LineStyle',':')
      end      

%% -------------- plot DECA evolution and pixel classification -------------- 
      figure(2);
      
      ycolor =[0 .5 1; .5 0 1; 0 1 1;1 0 1; 0 1 .5; 0 .5 .5; .5 1 0; 1 0 .5;.5 .5 .5; 1 1 .5];
      th=.99;
      beta = parameters_estimates.Post_probability;
      idx=find(prod(double(beta<=th)'));
      leg_cell = cell(1);
      plot(Y(1,idx),Y(2,idx),'.','Color',ycolor(1,:),'markersize',5)
      leg_cell{end} = 'data points';
      hold on
      for k = 1:k_opt,
          idx=find(beta(:,k)>th);
          if ~isempty(idx)
             plot(Y(1,idx),Y(2,idx),'.','Color',ycolor(1,:),'markersize',5);
             leg_cell{end+1} = 'data points';
        %     eval(['leg_cell{end +1} = ''class' int2str(1) ''';' ]);    
          end
      end
      plot(m_true(1,[1:p 1]), m_true(2,[1:p 1]),'k-','markersize',8,'Linewidth',2)
      leg_cell{end+1} = 'true ';
      plot(m_deca(1,[1:p 1]), m_deca(2,[1:p 1]),'g-', 'markersize',8,'markerfacecolor','y','Linewidth',2)
      leg_cell{end+1} = 'new algorithm';
      plot(squeeze(m_evolut(1,1:p,1)), squeeze(m_evolut(2,1:p,1)),'ro')
      leg_cell{end+1} = 'initial';
      
      plot(squeeze(m_evolut(1,1:p,:)), squeeze(m_evolut(2,1:p,:)),'y.')
      leg_cell{end+1} = 'evolution path';
      hold on;
      title('Data visualization') ;     
      legend(leg_cell);
     
%%/% -------------- plot DECA signatures -------------- 
      Md_aux = nan*zeros(length(wavlen),p);
      Md_aux = M_deca;
      
      p_color=['b','g','r','m'];
      figure(1)
      subplot(212)
      hold on
      for k=1:p
          % the results may need a permutation
          plot(wavlen,Md_aux(:,k),'--','Color',p_color(k));
      end
      axis([min(wavlen) max(wavlen) 0 1])
      hold off
      title('Endmembers Signatures and DECA estimates')
      xlabel('Wavelength (\mu m)')
      ylabel('Reflectance')
      
%% -
  
          
