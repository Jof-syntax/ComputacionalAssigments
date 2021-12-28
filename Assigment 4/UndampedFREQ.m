 function [MODES, FREQ] = UndampedFREQ(M,K,neig) 
if nargin == 0
    load('tmp.mat')
end

% We turn Mll and Kll symmetric, because otherwise the modes are not
%orthogonalized with respect to M. It should be noted that the lack
%of symmetry is caused by finine machine precision:   theoretically
%speaking, M and K are symmetric by construction
M = (M+M')/2; 
K = (K+K')/2;

  [MODES, EIGENVAL] = eigs(K,M,neig,'sm') ; 
  EIGENVAL = diag(EIGENVAL) ; 
  FREQ = sqrt(EIGENVAL) ;  
  
  [FREQ,imodes] = sort(FREQ) ;
  MODES = MODES(:,imodes); 
 end