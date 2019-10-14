function [MPost,Ctilde_all,W_coef]=fitGaussiansMixture(Xsp,fpost)
% 
% fit a mixture of multivariate Gaussians via BIC and EM_GM algorithm
% INPUT::
%  Xsp :: samples of multiple variables (n_sample \times n_variables)
% OUTPUT:: 
% Mpost :: the corresponding mean(s)
% Ctilte_all :: the covariance matrix(ces)
% W  :: corresponding weights
%
%  Note:: fpost is the minus log posterior fction handle for re-ordering the different
%  Gaussians (from most to min probable) 
%
tol = 1.e-2;

tol_keep_1=2.e-2;

k_s=length(Xsp(:,1));

% with at most 6 Gaussian mixtures
if length(Xsp(1,:)) < 6 
    k=length(Xsp(1,:));
else
    k=6;
end

stfd=0;  BIC_k_1=Inf;
while ~stfd && k>1 
    
    [W_a,M_a,V_a,L,E] = EM_GM(Xsp,k,[],[],[],[])  ;
    BIC_k=-2*L+k*log(length(k_s));
    fprintf('BIC for %i Gaussians : %f \n',k,BIC_k);
  
    if  BIC_k<BIC_k_1 || abs(BIC_k/BIC_k_1-1) < tol
        % ENSURE SUFFICIENT INCREASE of BIC to stop with previous estimate
        
        k=k-1;
        BIC_k_1=BIC_k;
        W=W_a;M=M_a; V=V_a; 
        
    else

        k=k+1;
        stfd=1; % BIC re-increases, keep the previous results

    end
    
end


if k==1          
    disp(['Only 1 Gaussian']);
    k=k+1; % reincrease to 
end
% sort the results in an ordered way, from the highest posterior value 
% this is where we particularize for function MinusLogPosterior etc.
m_log=zeros(k,1)
for p=1:k
        
        m_log(p)=feval(fpost,M(:,p));          
end

[m_log,IX]=sort(m_log,'ascend'); 
MPost=M(:,IX);
CC_all=V(:,:,IX);       % the associated fitted covariance matrix.
W_coef=W(IX);

if k==2  % case of 2 pdfs fitted
    %%% checks that the 2 pdf are actually sufficiently distinct
    
    dis=median(abs((M(:,1)-M(:,2))./mean(M')'));  % median here to relax bias
    
    if dis < tol_keep_1  % the 2 modes are really close to one another (2%) : JUST MERGE
        
        MPost=M(:,1);
        CC_all=V(:,:,1);
        W_coef=W(1)*0+1;
        k=1; % only keep ONE Gaussian in that case
    end
    
end

Ctilde_all=V;

end