function x=psi_inv(y)
%
% This function computes de inverse of 
% PSI function (digamma), i.e, y=psi(x)
% where x is real nonnegative.
% The function use the Newton's method
% with five iterations which are sufficient 
% to reach fourteen digits of precision.
%

% initial aproximation
% base on
% y=log(x-.5), if x>=.6
%  =-1./x-psi(1), if x<.6

    x = (y>=-2.22).*(exp(y)+.5) + ...;
        (y< -2.22).*(-1./(y-psi(1)));

for i=1:5
    x=x-(psi(0,x)-y)./psi(1,x);
end