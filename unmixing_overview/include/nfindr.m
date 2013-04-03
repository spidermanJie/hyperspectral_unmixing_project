function [Me, pp_indices, Yp] = nfindr_v1(Y,p)
%
% NFINDR Algorithm
%
% M. E. Winter, “N-FINDR: An algorithm for fast autonomous spectral
% endmember determination in hyperspectral data,” in Proc. SPIE Image
% Spectrometry V, 1999, vol. 3753, pp. 266–277.
%
% ------- Input variables -------------------------------------------
%
%  Y - matrix with dimensions L(channels) x N(pixels)
%      Each pixel is a linear mixture of p endmembers
%      signatures Y = M X, where M  and X are the mixing matrix
%      and the abundance fractions matrix, respectively.
%
%  p - number of endmembers in the scene
%
% ------- Output variables -------------------------------------------
%
% Me         - estimated mixing matrix (endmembers signatures)
%
% pp_indice - indexces of pixels chosen to be the most pure
%
% Yp         - Data Y projected on the identified signal subspace
%
%
% Author:  Jose M.  Bioucas Dias, January 2012


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(Y)
    error('there is no data');
else
    [L N]=size(Y);  % L number of bands (channels)
    % N number of pixels (LxC)
end

if (p<=0 | p>L | rem(p,1)~=0),
    error('ENDMEMBER parameter must be an  integer between 1 and L');
end


% orthogonal projection on an the best (p-1) affine set
y_m = mean(Y,2);
Corr = Y*Y'/N;
Cov  = Corr - y_m*y_m';
[U,S] = svds(Cov,p-1);         % computes the a (p-1)-orth basis
Y_o = Y - repmat(y_m,[1 N]);   % remove mean
Yp =  U' * Y_o;               % project the zero-mean data onto a 
                               % (p-1) subspace


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look for pixels in a random order
rand_index = randperm(N);
% pure pixels indeces
pp_indices = rand_index(1:p);

% increase dimension from p-1 to p
Yp = [ones(1,N); Yp];
M = Yp(:,pp_indices);


vol_simplex = 0;
for i=p+1:N
    y = Yp(:,rand_index(i));
    for j=1:p
        M_aux = M; 
        M_aux(:,j) = y;
        vol_aux(j) = abs(det(M_aux));
    end       
    [max_vol pos] = max(vol_aux);
    if vol_simplex < max_vol, 
        M(:,pos) = y; 
        vol_simplex = max_vol;
        pp_indices(pos) = rand_index(i);
    end
end



%  original data projected on the indentified subspace
Yp =  U * Yp(2:end,:) + repmat(y_m,[1 N]);   % again in dimension L


Me = Yp(:,pp_indices);  










