% data explore 

%endmembers
load USGS_pruned_10_deg; 

  [L n_materiais]=size(M);
   
    for  i=1:n_materiais 
   
    subplot(8,8,i), plot(M(:,i)); 
    hold on 
    end 

    % plot all endmembers in one figure
  parallelcoords(M')
  
  %real data(mixture) plot 
  
  
close all
clear all
verbose = 'off';

%% load data set  
load  Rcuprite

[B,n] = size(Y);

% estimate noise
[w Rw]=estNoise(Y);

% p == subspace dimension
% Ek == subspace orthogonal basis
[p,Ek] = hysime(Y,w,Rw,verbose);

% project onto the signal subspace
X = Ek'*Y; 


no_dims = round(intrinsic_dim(X', 'MLE'));
	disp(['MLE estimate of intrinsic dimensionality: ' num2str(no_dims)]);
	[mappedX, mapping] = compute_mapping(X', 'PCA', no_dims);	
	figure, scatter(mappedX(:,1), mappedX(:,2)); title('Result of PCA');
    gkde2(mappedX(:,1:2)');
    [mappedX, mapping] = compute_mapping(X, 'Laplacian', no_dims, 7);	
	figure, scatter(mappedX(:,1), mappedX(:,2), 5, labels(mapping.conn_comp)); title('Result of Laplacian Eigenmaps'); drawnow


