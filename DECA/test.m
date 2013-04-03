

clear all,
clc
close all
%% -------------- Simulation Options -------------- 
% option = 1 runs DECA on simulated data
% option = 2 runs DECA on Cuprite dataset

option = 1;

%% -------------- Simulation Parameters -------------- 

 
        p = 3;  
        % number of endmembers
        mm =3; %number of mixtures 
        N_LINES=50;                 % number of lines of the image
        N_COLUMNS=60;               % number of columns of the image 
        SNR = 2000;                  % signal-to-noise ratio (E ||y||^2/E ||n||^2) in dBs
        %SHAPE_PARAMETER = [4/10 6 25 9;
        %                   5/10 7 8 23
        %                   1/10 20 20  4];
        cluster = 2; 
%         weight_parameter = [3/5;2/5]; 
%         a = 0.5*ones(1,p); 
%         concentration = 5;
%         shape_pa =concentration*dirichlet(a,cluster); 
%         SHAPE_PARAMETER = [weight_parameter shape_pa] 
                          
        SHAPE_PARAMETER =[0.6    0.8    0.8    0.2;
                           0.4    0.4   0.4   0.5];
                                       
                           
        Kmod_max = cluster; 
        Kmod_min = cluster; 
                       
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
        clear datalib names aux st n_materiais sel_mat ;

       %% -------------- Generate the Hyperspectral Observations -------------- 
       % Observations are linear mixtures of the endmembers signatures plus noise

        N = N_LINES * N_COLUMNS;      % number of pixels
        [Y,x,noise,sigma,outliers,mode_length] = spectMixGen(M,N,'Source_pdf','Diri_mix','pdf_pars',SHAPE_PARAMETER,'snr',SNR);

        aux = cumsum([1; mode_length]);
        real_Kmod = length(mode_length);
        for i=1:real_Kmod
            eval(['Y' int2str(i) '=Y(:,aux(' int2str(i) '):aux(' int2str(i+1) ')-1);']);
        end % end of switch

%% --------------  Estimate noise matrix, assuming additive noise -------------- 
        [Rw w]=HyperNoise(Y,'additive','off'); 
        % Estimate number of endmembers
        HySime_Parameters = struct('Noise',w,'Noise_Correlation_Matrix',Rw, ...
                                   'Display_figure','off','Verbose','off');
        [p_hat,U,out_struct]=HySime(Y,HySime_Parameters);   
        U_ss = out_struct.Projection_Matrix_Percent(:,1:p);
        clear w Rw p_hat U

%% -------------- Project the dataset to the signal subspace -------------- 
        Y = U_ss'*Y;
        m_true =  U_ss'*M;
        % only for presentation purposes 
        if option==1
          for i=1:real_Kmod
              eval(['Y' int2str(i) '=U_ss''*Y' int2str(i) ';']);
          end
        end

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
      title(' prior distribution ');