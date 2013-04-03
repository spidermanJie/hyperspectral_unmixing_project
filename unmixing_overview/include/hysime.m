function [varargout]=hysime(varargin);
%
% HySime: Hyperspectral signal subspace estimation introduced in
%
% J. Bioucas-Dias and J. Nascimento, ``Hyperspectral subspace identification,'' 
% IEEE Transactions on Geoscience and Remote Sensing}, vol.~46., no.~8, 
% pp~ 2435-2445, 2005.
%
% -------------------------------------------------------------------------
%  Usage:
%
% [kf,Ek,E,delta_p]=hysime(y,w,Rw,verbose);
%
% -------------------------------------------------------------------------
% Input:
%        y -> hyperspectral data set (each column is a pixel)
%             with (L x N), where L is the number of bands
%             and N the number of pixels
%
%        w -> (L x N) matrix with the noise in each pixel Rw noise 
%              correlation matrix (L x L)
%
%        verbose ->  [optional] (on)/off
%
% -------------------------------------------------------------------------
%
% Output
%        kf -> signal subspace dimension
%
%        Ek -> matrix which columns are the eigenvectors that span 
%              the signal subspace
%
%        E -> Complete orthonormal basis for R^L ordered by 
%             delta_P = [Px(i)-Pn(i), i=1,...,L], where 
%             Px(i) and Pn(i) are, respectively, the noise and the signal
%             Power along the i-th eigenvector of Rx
%
%        delta_P -> delta_P = [Px(i)-Pn(i), i=1,...,L];
%
%
% NOTE: Ek = E(:,delta_P > 0).
%
% Authors: José Nascimento (zen@isel.pt)
%          José Bioucas-Dias (bioucas@lx.it.pt)
%
%
% -----------------------------------------------------------------------
% Copyright (2008): José Nascimento and Jose Bioucas-Dias
% 
% spectMixGen is distributed under the terms of 
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



error(nargchk(3, 4, nargin))
if nargout > 4 
    error('too many output parameters'); 
end
verbose = 1; % default value

y = varargin{1}; % 1st parameter is the data set
[L N] = size(y);
if ~numel(y)
    error('the data set is empty');
end
n = varargin{2}; % the 2nd parameter is the noise
[Ln Nn] = size(n);
Rn = varargin{3}; % the 3rd parameter is the noise correlation matrix
[d1 d2] = size(Rn);
if nargin == 4, 
    verbose = ~strcmp(lower(varargin{4}),'off');
end

if Ln~=L | Nn~=N,  % n is an empty matrix or with different size
   error('empty noise matrix or its size does not agree with size of y\n'),
end
if (d1~=d2 | d1~=L)
   fprintf('Bad noise correlation matrix\n'),
   Rn = n*n'/N; 
end    


x = y - n;

if verbose
    fprintf(1,'Computing the correlation matrices\n');
end
[L N]=size(y);
Ry = y*y'/N;   % sample correlation matrix 
Rx = x*x'/N;   % signal correlation matrix estimates 
if verbose
    fprintf(1,'Computing the eigen vectors of the signal correlation matrix\n');
end
[E,D]=svd(Rx); % eigen values of Rx in decreasing order, equation (15)
dx = diag(D);

if verbose
    fprintf(1,'Estimating the number of endmembers\n');
end
Rn=Rn+sum(diag(Rx))/L/10^5*eye(L);

Py = diag(E'*Ry*E); %equation (23)
Pn = diag(E'*Rn*E); %equation (24)
cost_F = -Py + 2 * Pn; %equation (22)
kf = sum(cost_F<0);
[cost_F_ascend,ind_asc] = sort( cost_F ,'ascend');
Ek = E(:,ind_asc(1:kf));
E = E(:,ind_asc);
if verbose
    fprintf(1,'The signal subspace dimension is: k = %d\n',kf);
end


if verbose
    indice=1:50;
    % only for plot purposes, equation (19)
    Py_sort =  trace(Ry) - cumsum(Py(ind_asc));
    Pn_sort = 2*cumsum(Pn(ind_asc));
    cost_F_sort = Py_sort + Pn_sort;
    
    figure(1000)
    set(gca,'FontSize',12,'FontName','times new roman')
    semilogy(indice,cost_F_sort(indice),'-',indice,Py_sort(indice),':', ...
        indice,Pn_sort(indice),'-.', 'Linewidth',2,'markersize',5)
    xlabel('k');
    ylabel('mse(k)');
    title('HySime')
    legend('Mean Squared Error','Projection Error','Noise Power')
end


varargout(1) = {kf};
if nargout >= 2
    varargout(2) = {Ek};
end
if nargout >= 3
    varargout(3) = {E};
end
if nargout >= 4
    varargout(4) = {-cost_F_ascend};
end

return
%end of function [varargout]=hysime(varargin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

