%% Name: geo_unmixing_npp
%
%   Generate the non-pure pixel unmixing results for section IV, Fig
%   15,top-right
%  
%   Data sets: 
%   SusgsP5SNR30 -> P = 5, non-pure pixels, SNR = 30 dB
%
%  Algorithms: VCA, NFINDR SISAL, MVC_NMF   
%  
%   [VCA]
%   J. M. P. Nascimento and J. M. Bioucas-Dias, “Vertex component
%   analysis: A fast algorithm to unmix hyperspectral data,” IEEE Trans.
%   Geosci. Remote Sens., vol. 43, no. 4, pp. 898–910, 2005.
%
%
%   [NFINDR] 
%   M. E. Winter, “N-FINDR: An algorithm for fast autonomous spectral
%   endmember determination in hyperspectral data,” in Proc. SPIE Image
%   Spectrometry V, 1999, vol. 3753, pp. 266–277
%
%
%   [SISAL]
%   J. Li and J. Bioucas-Dias, “Minimum volume simplex analysis: A
%   fast algorithm to unmix hyperspectral data,” in Proc. IEEE Int. Conf.
%   Geosci. Remote Sens. (IGARSS), 2008, vol. 3, pp. 250–253.
%
%
%   [MVC_NMF]
%   L. Miao and H. Qi, "Endmember extraction from highly mixed data
%   using minimum volume constrained nonnegative matrix factorization,"
%   IEEE Transactions on Geoscience and Remote Sensing, vol. 45, 
%   no. 3, pp. 765–777, 2007.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MVC_NMF is not  included  in this version
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Jose M. Bioucas-Dias (bioucas@lx.it.pt), January, 2012

%% begining script

clear all,
%clc
close all

randn('seed',200)
rand('seed',200)


% ----- pure pixels ----------------
load SusgsP5SNR30

[B,n] = size(Y);
% number of endmembers
p = 5;

% Compute SNR-SD to evaluate the dificulty of the SU problem
% estimate noise
[w,Rw] = estNoise(Y);

% estimate signal subspace
[phat,E] = hysime(Y,w,Rw,'off');

Rx = Y*Y'/n - Rw;
SNR_SD = diag(E'*Rx*E)./diag(E'*Rw*E);

% compute in dB
SNR_SD = 10*log10(SNR_SD(1:p))

fprintf('\nMinimum SNR_SD %3.2f \n',min(SNR_SD))

% Project data on the signal subspace
Y = E*E'*Y;


%%
%--------------------------------------------------------------------------
%       Project  on the affine set defined by the data in the l2 sense
%-------------------------------------------------------------------------
%
%   The application of this projection ensures that the data belongs to
%   an affine set.
%
%   Up is an isometric matrix that spans the subspace where Y lives
%   Yp contains the coordinates wrt Up

[Yp,Up,my,angles,scales] = affineProj(Y,p,'proj_type','orth');



%%
%--------------------------------------------------------------------------
%         ALGORITHMS
%-------------------------------------------------------------------------


% set which algorithms will run
run_vca    = 1;
run_nfindr = 1;
run_sisal  = 1;
run_mvc    = 0;
run_spice  = 0;
run_deca   = 0;



%%
%--------------------------------------------------------------------------
%         VCA  - Vertex component analysis
%-------------------------------------------------------------------------
%
% start timer

t(1) = inf;
if run_vca == 1
    % start timer
    tic
    [Mvca, indice,ys]= VCA(Yp,'Endmembers',p,'SNR',0);
    % stop timer
    t(1) = toc;
end
%%
%--------------------------------------------------------------------------
%         NFINDR -
%-------------------------------------------------------------------------
%

t(2) = inf;
if run_nfindr == 1
    % start timer
    tic
    [Mnf, indice,dummy]= nfindr(Yp,p);
    % stop timer
    t(2) = toc;
end
%%
%--------------------------------------------------------------------------
%         SISAL[1] -  Simplex identification via split augmented Lagrangian
%-------------------------------------------------------------------------

t(3) = inf;
if run_sisal == 1
    % start timer
    tic
    % set the hinge regularization parameter
    tau = 0.03;
    [Msisal] = sisal(Yp,p, 'spherize', 'no','MM_ITERS',80, 'TAU',tau, 'verbose',2);
    drawnow;
    t(3)=toc;
end



%%
%--------------------------------------------------------------------------
%         MVC-NMF 
%-------------------------------------------------------------------------

t(4)=inf;
if run_mvc  == 1
    lambda = 0.0001;
    % SUNSAL
    % initial abundances determined with the Mvca mining matrix 
    [sinit] = sunsal(Mvca,Yp,'lambda',lambda,'ADDONE','yes','POSITIVITY','yes', ...
        'AL_iters',2000,'verbose','yes');

    [PrinComp, meanData] = pca((Up*Yp)', 0);

    % MVCNMF
    tol = 1e-6;
    maxiter = 150;
    T = 0.015;
    showflag = 1;
    tic
    % use conjugate gradient to find A can speed up the learning
    [Mmvc, sest] = mvcnmf(Up*Yp,Up*Mvca,sinit,M,Up,PrinComp,meanData,T,tol,maxiter,showflag,2,1);
    Mmvc = Up'*Mmvc;
    t(4)=toc;
end


%%
%--------------------------------------------------------------------------
%         DECA
%-------------------------------------------------------------------------

t(5) = inf;
if run_deca == 1
tic;
    
    Kmod_max=1;
    Kmod_min=1;
    max_iter=300;
    th=1e-5;
    ver_fig_deca = 1;
    verbose = 1;

    [M_deca,k_opt,a_opt,e_opt,L_global,e_global,transitions]= ...
        deca_mdl(Yp,p,Kmod_max,Kmod_min,max_iter,th,Up,ver_fig_deca,verbose);

t(5) = toc;   
end


%%
%--------------------------------------------------------------------------
%         SPICE
%-------------------------------------------------------------------------
t(6) = inf;
if run_spice == 1
    tic;
    parameters.u = 0.001; %Trade-off parameter between RSS and V term
    parameters.gamma = 0.5; %Sparsity parameter
    parameters.M = p; %Initial number of endmembers
    parameters.endmemberPruneThreshold = -inf;  % this means that NO PRUNING
    parameters.changeThresh = 1e-4; %Used as the stopping criterion
    parameters.iterationCap = 300; %Alternate stopping criterion
    parameters.produceDisplay = 1; %verbose;
    parameters.initEM =  Mvca; %nan; %This randomly selects parameters.M initial endmembers from the input data
    [Mspice,Sspice] = SPICE(Yp, parameters);

    t(6) = toc;
end



%%
%--------------------------------------------------------------------------
%         Project the original mixing matxix and the data set the
%         identified affine set.
%-------------------------------------------------------------------------
MT = M;
Mtrue = Up'*M;
Y=Up'*Y;


%%
%--------------------------------------------------------------------------
%        Display the results
%-------------------------------------------------------------------------

% selects axes  to display
%

I = 1;
J = 2;
K = 3;

% canonical orthogonal directions
E_I = eye(p);

v1 = E_I(:,I);
v2 = E_I(:,J);
v3 = E_I(:,K);

% original axes

Q = inv(Mtrue);
% v1 = Q(I,:)';
% v2 = Q(J,:)';
% v3 = Q(K,:)';

Y = [v1 v2 v3]'*Y;
m_true = [v1 v2 v3]'*Mtrue;



% legend
leg_cell = cell(1);
leg_cell{end} = 'data points';
H_2=figure;
plot(Y(1,:),Y(2,:),'k.','Color',[ 0.2 0.2 0.2])

hold on;
plot(m_true(1,[1:p 1]), m_true(2,[1:p 1]),'o', 'Color',[0 0 0],'LineWidth',2)
leg_cell{end +1} = 'true';

if run_nfindr
    m_findr  = [v1 v2 v3]'*Mnf;
    plot(m_findr(1,[1:p 1]), m_findr(2,[1:p 1]),'v', 'Color',[0  0 1],'LineWidth',2)
    leg_cell{end +1} = 'NFINDR';
end


if run_vca == 1
    m_vca  = [v1 v2 v3]'*Mvca;
    plot(m_vca(1,[1:p 1]), m_vca(2,[1:p 1]),'^', 'Color',[0  0.8 0],'LineWidth',2)
    leg_cell{end +1} = 'VCA';
end

if run_mvc
    m_mvc  = [v1 v2 v3]'*Mmvc;
    plot(m_mvc(1,[1:p 1]), m_mvc(2,[1:p 1]),'^', 'Color',[1  0.5 0],'LineWidth',2)
    leg_cell{end +1} = 'MVC-NMF';
end


if run_sisal == 1
    m_sisal = [v1 v2 v3]'*Msisal;
    plot(m_sisal(1,1:p), m_sisal(2,1:p),'S', 'Color',[1 0.0 0],'LineWidth',2)
    leg_cell{end +1} = 'SISAL';
end
  

if run_spice == 1
     m_spice  = [v1 v2 v3]'*Mspice;
     plot(m_spice(1,[1:p 1]), m_spice(2,[1:p 1]),'v', 'Color',[0.5  0.5 0.5],'LineWidth',2)
     leg_cell{end +1} = 'SPICE';
 end

 if run_deca == 1
     m_deca  = [v1 v2 v3]'*M_deca;
     plot(m_deca(1,[1:p 1]), m_deca(2,[1:p 1]),'>', 'Color',[0.7  0.2 0.2],'LineWidth',2)
     leg_cell{end +1} = 'SPICE';
 end

    
xlabel('e_1^T Y^T'),ylabel('e_2^T Y^T');
legend(leg_cell)
title('Endmembers and data points (2D projection) - pure pixel','FontSize',12)

axis off
set(gca,'FontSize',12)



fprintf('\nTIMES (sec):\n VCA = %3.2f\n NFINDR = %3.2f\n SISAL = %3.2f\n MVC = %f\n SPICE = %3.2f\n DECA = %3.2f\n\n', ...
     t(1), t(2),t(3),t(4),t(5),t(6))

%% 
%--------------------------------------------------------------------------
%        Display errors
%-------------------------------------------------------------------------
% alignament

VCA_ERR = inf;
VCA_MAX_ANG = inf;
if run_vca
    angles = Mtrue'*Mvca./(repmat(sqrt(sum(Mtrue.^2)),p,1)'.*(repmat(sqrt(sum(Mvca.^2)),p,1)));
    P = zeros(p);
    for i=1:p
        [dummy,j] = max(angles(i,:));
        P(j,i) = 1;
        angles(:,j) = -inf;
    end
    % permute colums
    Mvca = Mvca*P;

    VCA_ERR =norm(Mtrue-Mvca,'fro')/norm(Mtrue,'fro');
    % maximum angle
    angles = sum(Mtrue.*Mvca)./sqrt(sum(Mtrue.^2).* sum(Mvca.^2));
    VCA_MAX_ANG = acos(min(angles(:)))*180/pi;
end


FINDR_ERR = inf;
FINDR_MAX_ANG = inf;
if run_nfindr
    angles = Mtrue'*Mnf./(repmat(sqrt(sum(Mtrue.^2)),p,1)'.*(repmat(sqrt(sum(Mnf.^2)),p,1)));
    P = zeros(p);
    for i=1:p
        [dummy,j] = max(angles(i,:));
        P(j,i) = 1;
        angles(:,j) = -inf;
    end
    % permute colums
    Mnf = Mnf*P;
    FINDR_ERR =norm(Mtrue-Mnf,'fro')/norm(Mtrue,'fro');
    
    % maximum angle
    angles = sum(Mtrue.*Mnf)./sqrt(sum(Mtrue.^2).* sum(Mnf.^2));
    FINDR_MAX_ANG = acos(min(angles(:)))*180/pi;
end




SISAL_ERR = inf;
SISAL_MAX_ANG = inf;
if run_sisal
    angles = Mtrue'*Msisal./(repmat(sqrt(sum(Mtrue.^2)),p,1)'.*(repmat(sqrt(sum(Msisal.^2)),p,1)));
    P = zeros(p);
    for i=1:p
        [dummy,j] = max(angles(i,:));
        P(j,i) = 1;
        angles(:,j) = -inf;
    end
    % permute colums
    Msisal = Msisal*P;
    SISAL_ERR =norm(Mtrue-Msisal,'fro')/norm(Mtrue,'fro');
   
    % compute angles 
    angles = sum(Mtrue.*Msisal)./sqrt(sum(Mtrue.^2).* sum(Msisal.^2));
    SISAL_MAX_ANG = acos(min(angles(:)))*180/pi;
end


MVC_ERR = inf;
MVC_MAX_ANG = inf;
if run_mvc == 1
    angles = Mtrue'*Mmvc./(repmat(sqrt(sum(Mtrue.^2)),p,1)'.*(repmat(sqrt(sum(Mmvc.^2)),p,1)));
        MVC_MAX_ANG = acos(min(angles(:)))*180/pi;
    P = zeros(p);
    for i=1:p
        [dummy,j] = max(angles(i,:));
        P(j,i) = 1;
        angles(:,j) = -inf;
    end
    % permute colums
    Mmvc = Mmvc*P;

    MVC_ERR =norm(Mtrue-Mmvc,'fro')/norm(Mtrue,'fro');
    
    % compute angles 
    angles = sum(Mtrue.*Mmvc)./sqrt(sum(Mtrue.^2).* sum(Mmvc.^2));
    MVC_MAX_ANG = acos(min(angles(:)))*180/pi;
end

SPICE_ERR = inf;
SPICE_MAX_ANG = inf;
if run_spice == 1
    angles = Mtrue'*Mspice./(repmat(sqrt(sum(Mtrue.^2)),p,1)'.*(repmat(sqrt(sum(Mspice.^2)),p,1)));
    P = zeros(p);
    for i=1:p
        [dummy,j] = max(angles(i,:));
        P(j,i) = 1;
        angles(:,j) = -inf;
    end
    % permute colums
    Mspice = Mspice*P;

    SPICE_ERR =norm(Mtrue-Mspice,'fro')/norm(Mtrue,'fro');
    
    % compute angles 
    angles = sum(Mtrue.*Mspice)./sqrt(sum(Mtrue.^2).* sum(Mspice.^2));
    SPICE_MAX_ANG = acos(min(angles(:)))*180/pi;
end

DECA_ERR = inf;
DECA_MAX_ANG = inf;
if run_deca == 1
    angles = Mtrue'*M_deca./(repmat(sqrt(sum(Mtrue.^2)),p,1)'.*(repmat(sqrt(sum(M_deca.^2)),p,1)));

    P = zeros(p);
    for i=1:p
        [dummy,j] = max(angles(i,:));
        P(j,i) = 1;
        angles(:,j) = -inf;
    end
    % permute colums
    M_deca = M_deca*P;

    DECA_ERR =norm(Mtrue-M_deca,'fro')/norm(Mtrue,'fro');
    %compute angles 
    angles = sum(Mtrue.*M_deca)./sqrt(sum(Mtrue.^2).* sum(M_deca.^2));
    DECA_MAX_ANG = acos(min(angles(:)))*180/pi;
end



st = sprintf('ERROR(mse):\n VCA = %f\n NFINDR = %f\n SISAL = %f\n MVC-NMF = %f\n SPICE = %f\n DECA = %f\n', ...
   VCA_ERR, FINDR_ERR,  SISAL_ERR,  MVC_ERR, SPICE_ERR, DECA_ERR);
fprintf(strcat('\n ', st, '\n'));



st = sprintf('MAX_ERROR(deg):\n VCA = %f\n NFINDR = %f\n SISAL = %f\n MVC-NMF = %f\n SPICE = %f\n DECA = %f\n', ...
   VCA_MAX_ANG, FINDR_MAX_ANG,  SISAL_MAX_ANG,  MVC_MAX_ANG, SPICE_MAX_ANG, DECA_MAX_ANG);
fprintf(strcat('\n ', st, '\n'));



%--------------------------------------------------------------------------
%        Plot signatures
%-------------------------------------------------------------------------
% Choose signatures

leg_cell = cell(1);
H_3=figure;
hold on
clear p_H;

% plot signatures
p_H(1) = plot(1:B,(MT(:,1))','k','LineWidth',2);
leg_cell{end} = 'Mtrue';
if run_sisal == 1
    p_H(2) =plot(1:B,(Up*Msisal(:,1))','r','LineWidth',2);
    leg_cell{end+1} = 'Msisal';
end
if run_vca == 1
    p_H(3) =plot(1:B,(Up*Mvca(:,1))','g','LineWidth',2);
    leg_cell{end+1} = 'Mvca';
end
if run_nfindr
    p_H(4) =plot(1:B,(Up*Mnf(:,1))','b','LineWidth',2);
    leg_cell{end+1} = 'Mnfindr';
end

if run_mvc == 1
p_H(5) =plot(1:B,(Up*Mmvc(:,1))','m','LineWidth',2);
leg_cell{end+1} = 'Mmvc-nmf';
end

if run_spice == 1
    p_H(6) =plot(1:B,(Up*Mspice(:,1))','c','LineWidth',2);
    leg_cell{end+1} = 'Mspice';
end

if run_deca == 1
    p_H(7) =plot(1:B,(Up*M_deca(:,1))','k.','LineWidth',2);
    leg_cell{end+1} = 'Mdeca';
end


for i=2:p
    plot(1:B,(MT(:,i))','k','LineWidth',2);
    if run_sisal == 1
        plot(1:B,(Up*Msisal(:,i))','r','LineWidth',2);
    end
    if run_vca == 1
        plot(1:B,(Up*Mvca(:,i))','g','LineWidth',2);
    end
    if run_nfindr == 1
        plot(1:B,(Up*Mnf(:,i))','b','LineWidth',2);
    end
    if run_mvc == 1
        plot(1:B,(Up*Mmvc(:,i))','m','LineWidth',2);
    end
    if run_spice == 1
        plot(1:B,(Up*Mspice(:,i))','c','LineWidth',2);
    end
    if run_deca == 1
        plot(1:B,(Up*M_deca(:,i))','k.','LineWidth',2);
    end
end

legend(leg_cell)
title('Endmembers')
xlabel('spectral band')

pos2 = get(H_2,'Position');
pos2(1)=50;
pos2(2)=1;
set(H_2,'Position', pos2)

pos3 = get(H_3,'Position');
pos3(1)=600;
pos3(2)=100+400;
set(H_3,'Position', pos3)



%--------------------------------------------------------------------------
%        Plot the 1st signatures
%-------------------------------------------------------------------------
% Choose signatures

leg_cell = cell(1);
H_3=figure;
hold on
clear p_H;

% plot signatures
p_H(1) = plot(1:B,(MT(:,1))','k','LineWidth',2);
leg_cell{end} = 'true';
if run_vca == 1
    p_H(3) =plot(1:B,(Up*Mvca(:,1))','g','LineWidth',2);
    leg_cell{end+1} = strcat('vca',sprintf(' (%1.1f)', VCA_MAX_ANG ));
end
if run_nfindr
    p_H(4) =plot(1:B,(Up*Mnf(:,1))','b','LineWidth',2);
    leg_cell{end+1} = strcat('nfindr',sprintf(' (%1.1f)', FINDR_MAX_ANG ));
end
if run_sisal == 1
    p_H(2) =plot(1:B,(Up*Msisal(:,1))','r','LineWidth',2);
    leg_cell{end+1} =  strcat('sisal',sprintf(' (%1.1f)', SISAL_MAX_ANG ));
end
if run_mvc == 1
p_H(5) =plot(1:B,(Up*Mmvc(:,1))','m','LineWidth',2);
leg_cell{end+1} = strcat('mvc-nmf',sprintf(' (%1.1f)', MVC_MAX_ANG ));
end

if run_spice == 1
    p_H(6) =plot(1:B,(Up*Mspice(:,1))','c','LineWidth',2);
    leg_cell{end+1} = strcat('spice',sprintf(' (%1.1f)', SPICE_MAX_ANG ));
end

if run_deca == 1
    p_H(7) =plot(1:B,(Up*M_deca(:,1))','k.','LineWidth',2);
    leg_cell{end+1} = trcat('deca',sprintf(' (%1.1f)', SPICE_MAX_ANG ));
end

% for i=2:p
%     plot(1:B,(MT(:,i))','k','LineWidth',2);
%     if run_sisal == 1
%         plot(1:B,(Up*Msisal(:,i))','r','LineWidth',2);
%     end
%     if run_vca == 1
%         plot(1:B,(Up*Mvca(:,i))','g','LineWidth',2);
%     end
%     if run_nfindr == 1
%         plot(1:B,(Up*Mnf(:,i))','b','LineWidth',2);
%     end
%     if run_mvc == 1
%         plot(1:B,(Up*Mmvc(:,i))','m','LineWidth',2);
%     end
%     if run_spice == 1
%         plot(1:B,(Up*Mspice(:,i))','c','LineWidth',2);
%     end
%     if run_deca == 1
%         plot(1:B,(Up*M_deca(:,i))','k.','LineWidth',2);
%     end
%     if run_sr == 1
%         plot(1:B,(Msr(:,i))','g.','LineWidth',2);
%     end
% end

legend(leg_cell)
title('1st spectral signature (error - deg)', 'Fontsize', 12)
xlabel('spectral band')

pos2 = get(H_2,'Position');
pos2(1)=50;
pos2(2)=1;
set(H_2,'Position', pos2)

pos3 = get(H_3,'Position');
pos3(1)=600;
pos3(2)=100+400;
set(H_3,'Position', pos3)

axis off

set(gca,'FontSize',12)



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %