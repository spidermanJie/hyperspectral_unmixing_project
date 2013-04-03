function   [Teps, Talpha, Source] = generateSamples(m,p,N); 

%Data = A*S 
%   TA   \\true A
%   TW    \\true W TW=inv(TA); 
%   Teps, Talpha




% 
% m is number of mixtures in the sources
% p is number of sources
% Source is N by p!!! not p by N. where N is the number of samples
% W is p by p W = inv(A); 
% eps is 1 by m ; 
% alpha is m by p; 




%p = 4; %number of sources
%m = 3; %number of mixtures
%N = 1000; % numnber of Samples
Talpha = randi([2,20],m,p) ; % the True parameter of Dir distribution
Teps = randi([2,20],1,m);% the True parameter for the mixture weights 

Teps = Teps/sum(Teps); 

%generate samples
r = mnrnd(1,Teps,N) ;
Source = zeros(N,p);
sum2 = sum(r); 

for i =1:m
    ind = find(r(:,i)); 
    S1 = sample_dirichlet(Talpha(i,:),sum2(i)); 
    Source(ind,:) = S1; 
    
end
 
Source =formatData(Source);