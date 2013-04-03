%%  The code and data herein distributed reproduce the results published in
%  the paper 
%
%  J.  Bioucas-Dias, A. Plaza, N. Dobigeon, M. Parente, Q. Du, P. Gader,  
%  and J. Chanussot, "Hyperspectral unmixing overview: geometrical, 
%  statistical, and sparse regression-based approaches",  IEEE Journal of 
%  Selected Topics in Applied Earth Observations and Remote Sensing, vol. 5, 
%  no. 2, pp. 354-379, 2012.
%  
%
%% The code is distributed in the folder "paper_scripts" and is organized 
%  according to the paper sections. Below is listed the correspondence 
%  between directories and sections and, for each directory,  the list of 
%  scripts it contains and the respective figures in the paper.
%  
%       paper_scripts\linear_mixing_model (Section II)
%
%            snr_sd.m ->  Fig. 8
%          
%       paper_scripts\subspace_identification (Section III)
%
%            projSpectra.m ->  Fig. 9
%
%            denoise_eigen_image.m ->  plots in Fig. 10 (See Note 1 below)
%
%            orth_dpft_proj.m ->   Fig. 12
%
%
%       paper_scripts\geometric_unmixing (Section IV)
%
%           geo_unmixing_pp.m    -> Fig. 15, top-left
%           geo_unmixing_npp.m   -> Fig. 15, top-right
%           geo_unmixing_mp08.m  -> Fig. 15, bottom-left
%           geo_unmixing_hmX10.m -> Fig. 15, bottom-right
%
%           NOTE: MVC-NMF  results are not included in this version
%
%
%      paper_scripts\statistical_methods
%
%          Fig. 16  ->  (TO BE INLUDED SOON)   
%          Fig. 17  ->  (TO BE INLUDED SOON)   
%
%      paper_scripts\sparse_regression (Section VI)
%
%           sparse_regression.m    -> Fig. 19
%
%            
%
%%  Notes:
%
%   1) Package instalation: unzip the files to a directory and run the setup
%      in this directory. The scripts listed above are then ready to run.
%
%
%   2) The script  denoise_eigen_image.m uses the package BM3D 
%      (v1.9 (26 August 2011))to implement the denoising algorithm BM3D 
%      introduced in 
%
%      K. Dabov, A. Foi, V. Katkovnik, and K. Egiazarian, "Image denoising by
%      sparse 3D transform-domain collaborative filtering," IEEE Trans. 
%      Image Process., vol. 16, no. 8, pp. 2080-2095, August 2007.
%
%      The BM3D package  is available at the 
%      BM3D web page:  http://www.cs.tut.fi/~foi/GCF-BM3D
%
%      Download this package and install it is the folder include/BM3D
%      
%
%   
%% ACKNOWLEDGMENTS
%
% The authors acknowledge the following individuals and organisations:
%
%   - Dr. L. Miao for providing the matlab code fot the MVC-NMF algorithm
%
%       L. Miao and H. Qi, "Endmember extraction from highly mixed data
%       using minimum volume constrained nonnegative matrix factorization,"
%       IEEE Transactions on Geoscience and Remote Sensing, vol. 45, 
%       no. 3, pp. 765–777, 2007.
%
%   - Dr. Jose Nascimento for his support on the DECA code and data.
%
%   - Dr. Robert O. Green and the AVIRIS team for making 
%     the Rcuprite hyperspectral data set available to the community
%
%   - The United States Geological Survey (USGS) for their publicly available 
%     library of mineral signatures. 
%
%   - The  Army Geospatial Center, US Army Corps of Engineers, 
%     for making the HYDICE  Terrain data set available to the community.
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Jose M. Bioucas Dias, March 2012
%

