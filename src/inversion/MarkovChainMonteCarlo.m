function [accept, Xf, PXf, cf ,k_o, stab]=MarkovChainMonteCarlo(fid_log, mstart,Cstart,Smax,plot_flag,AutomaticStop,fun,varargin)
% reversible jump
%   Metropolis Hastings Algorithm from Gelmann 1995
%
%  Author: Brice Lecampion
%  Project : Initial stress inversion, pdd # 61922
%  last modification : August 2009
%
%   INPUT VARIABLES
%  fid_log :: the fid of log file
%  mstart :: starting point of the chain
%  Cstart :: Covariance at the starting point of the chain
%  Smax :: max number of steps of the chain
%  plot_flag :: flag for plotting the pdf current val during the chain
%  progression (1) or not (0)
%  AutomaticStop :: flag for AutomaticStop of the chain (1) or not (0)
%  fun :: string function name with additional variables in varargin
%           COMPUTES -log(pdf)  (not pdf), calling sequence is fun(m, varargin)
%
%
%  OUTPUT
% accept : vector with acceptance for each attempted jumps (0 or 1)
% Xf : matrix containing parameters chain results
% PXf : vector of accepted -log posterior value
% cf : jump tuning factor (value at last jump)
% k_o : integer (define end of the burn in period: ie. when chain appears stabilized)
% stab :: integer flag (0, or 1) stabilized chain (1) or not (0)
%

if(fid_log==-1)
    error('invalid fid of log file');
end

q=length(mstart); %% dimension of the problem
L=chol(Cstart,'lower');
cf=2.4;
k=1;   % iterator for the acceptance
jad=1; % iterator for accept rate tuning
i=1;   % iterator of the chain

mt=mstart;
m_log_piT=feval(fun,mstart,varargin{:}); %% -log of likelihood
PX=zeros(Smax,1);
X=zeros(Smax+1, length(mstart));
accept=zeros(Smax,1);
X(1,:)=[mstart'];


%   stabilization criterium tolerance
%eps_s = 1e-4;		% initial setup used by TLei
eps_s = 1e-3;		% modified by Brice
%eps_s = 0.02;		% increased further

stab=0;
fprintf('Chain will be assumed stationary if the   mode with the max post pdf do not change by more than %i percent (in relative terms)\n', 100*eps_s);
fprintf(fid_log,'Chain will be assumed stationary if the mode with the max post pdf do not change by more than %i percent (in relative terms)\n', 100*eps_s);

%em stab
em_ng=4;
M_k_1=[X(1,:)'  X(1,:)' X(1,:)' X(1,:)' ];
W_k_1=[1 0 0 0];

%L_stab=4000*q;   % length of the chain when stationary state has been reached (i.e for sampling purpose)
L_stab = 5000*q;


PX_k_1=1;
warning off all;

%while (i<Smax) ;
while ( i <= Smax )

  if ( mod(i,1e4) == 0 )
    fprintf(1, '\n');
    fprintf(1, 'MCMC Step #%d,000\n', i/1000);
  end
    
    % JUMP
    
    w=randn(q,1);            %random('normal',0,1,q,1);
    Y=mt+(cf/sqrt(q))*(L*w);
    
    while ~all(isfinite(Y))
        
        w=randn(q,1);            %random('normal',0,1,q,1);
        Y=mt+(cf/sqrt(q))*(L*w);
        
    end
    
    % evaluate  - log of posterior
    m_log_piY=feval(fun,Y,varargin{:});
    r=exp(m_log_piT-m_log_piY); %%% siwtch back to proper ratio
    
    % ACCEPTANCE CRITERIA
    if ( r > 1)
        % accept
        
        mt=Y;
        m_log_piT=m_log_piY;
        accept(i)=1;
        
        X(k+1,:) =[mt'];
        PX(k)=m_log_piT;
        k=k+1;
        %  disp('accept ');disp(k);disp((mod(k,10*q)==0) && (accept(i)==1));
    else
        % conditional acceptance
        alpha=rand(1);%random('unif',0,1);
        if (alpha < r ) ;
            % accept
            
            mt=Y;
            m_log_piT=m_log_piY;
            accept(i)=1;
            X(k+1,:) =[mt'];
            PX(k)=m_log_piT;
            k=k+1; % increase accepted steps iterator
            % disp('accept c');disp(k);disp((mod(k,10*q)==0) && (accept(i)==1));
            
        else
            % reject
            accept(i)=0;
        end
    end
    
    % TUNE ACCEPTANCE RATE :: acceptance rate should remain between 0.2 and
    % 0.44
    if ( (jad/200) > 1)
        if (mean(accept(1:i))<0.2)
            cf=cf/2.;
        else
            if (mean(accept(1:i))>0.44)
                cf=cf*2.;
            end
        end
        jad=1; % re-initialization
    end
    
    
    % Draw ?
    if ( (mod(k,1000)==0) && (plot_flag==1))  
        plot(PX(1:k-1));drawnow;
    end
    
    %%%% Check stabilization of the chain
    
    if ( (mod(k,100*q)==0) && (accept(i)==1) ) % every 50*q accepted samples
        
        
        %%%% perform EM_algo on the samples to fit Gaussian mixtures of em_ng gaussian
        [W_k,M_k] = EM_GM(X(10:3:k-1,:),em_ng,[],[],[],[])  ;  % throw the first tenth seems to stabilize results
        
        % sort according to the first coordinates of the model parameters
        % vector
        
        [W_k,ix]=sort(W_k,'descend');
        
        M_k=M_k(:,ix);
        
        %%% compute the log posterior for all estimates, and make the test
        %%% on the lowest one.
        m_log=zeros(em_ng,1);
        for h=1:em_ng
            if all(isfinite(M_k(:,h)))
                m_log(h)=feval(fun,M_k(:,h),varargin{:});
            else
                m_log(h)=NaN;
            end
        end
        
        PX_k=min(m_log);    % always compare to the modes with the minimum posterior
        
        crit_p=norm(1-PX_k/PX_k_1);
        
        % find the previous estimated modes the closest to the two highest
        % (with the most weight)
        % current modes
        
%         [c,i1]=min([norm(M_k(:,1)-M_k_1(:,1)),norm(M_k(:,1)-M_k_1(:,2)),norm(M_k(:,1)-M_k_1(:,3)),norm(M_k(:,1)-M_k_1(:,4))]);
%         
%         crit_max_modes=norm(1-M_k(:,1)./M_k_1(:,i1));%+norm(1-M_k(:,2)./M_k_1(:,i2));
%         crit_Wk=norm(W_k-W_k_1);
        
        W_k_1=W_k;
        M_k_1=M_k;
        PX_k_1=PX_k;
        
        %   disp([' crit on X ',num2str(crit)]);
        %disp([' crit on PX ',num2str(crit_p)]);
        %disp([' crit on W_k ',num2str(crit_Wk)]);
        %disp([' crit on Max modes ',num2str(crit_max_modes) ] );
        %   disp([' crit on 3 highest modes ',num2str(crit_modes) ] );
        
        if (  crit_p < eps_s )
            
	    % first accepted
            if (stab==0)
                fprintf('Chain appears stabilized w/r to "mode" after  %i accepted steps,  criteria is = %f  percent, current mlogP =%f\n',k-1,crit_p*100,m_log_piT);
                fprintf(fid_log, 'Chain appears stabilized w/r to "mode" after  %i accepted steps, criteria is = %f percent, current mlogP =%f\n',k-1,crit_p*100,m_log_piT);
                stab=1;
                k_o=k-1; % iterator w
                
            else
                fprintf('Chain stationary  criteria after %i accepted steps is = %f percent, current mlogP =%f\n',k-1,crit_p*100,m_log_piT);
                fprintf(fid_log, 'Chain stationary  criteria %i accepted steps is = %f percent, current mlogP =%f\n',k-1,crit_p*100,m_log_piT);
            end
            
        else
            
            if (stab==1)  % CHAIN IS NO MORE STABILIZED, print and
                
                fprintf('Current chain is no more stationary  ! \n');
                fprintf(fid_log, 'Current chain is no more stationary ! \n');
                fprintf('Chain stationary  criteria after %i accepted steps is = %f percent, current mlogP =%f\n',k-1,crit_p*100,m_log_piT);
                fprintf(fid_log, 'Chain stationary  criteria after %i accepted steps is = %f percent, current mlogP =%f\n',k-1,crit_p*100,m_log_piT);
                stab=0;
                
            else
                
                %fprintf('Chain stationary  criteria after %i accepted steps is  = %f percent, current mlogP =%f\n',k-1,crit_p*100,m_log_piT);
                fprintf('Chain stationary  criteria after %i accepted steps is = %f percent, current mlogP =%f\n',k-1,crit_p*100,m_log_piT);
                fprintf(fid_log, 'Chain stationary  criteria after %i accepted steps is = %f percent, current mlogP =%f\n',k-1,crit_p*100,m_log_piT);
                
            end
            
        end
        
    end
    
    % case chain is stabilized , test for the stoppting critera (i.e.
    % chain stabilized and has run for L_stab step in the stabilized
    % state
    
    if  ( (stab==1)  && ((k-1-k_o)>L_stab) && (AutomaticStop==1) )
                 fprintf('End of MCMC : chain appears stationary and   sampled ');
                 fprintf(fid_log, 'End of MCMC : chain appears stationary and sampled ');
        break
    end
    
    
    % iterators increase
    jad=jad+1;
    i=i+1;
    
    
end
%end while

%%% PRINT DIAGNOSTICS
if (stab ~=  1 )
    
    fprintf('Chain has not reached a stationary state after the specified number of steps: Dubious results, run the chain longer by oversetting default ');
    fprintf(fid_log, 'Chain has not a stationary state after the specified number of steps: Dubious results, run the chain longer by oversetting default');
    k_o=k-1;
    
else
    
    % stationary chain
    if ((k-1-k_o)>L_stab) % has been fully sampled
        
        fprintf('Chain appears stationary and sampled after %i steps of which %i have been accepted',i,k);
        fprintf(fid_log,'Chain appears stationary and sampled after %i steps of which %i have been accepted',i,k);
        
    elseif ( (i==Smax)  ) && ( (k-1-k_o)<L_stab  )  % max number of steps and not "fully" sampled        % check for sufficient sampling ......even though chain appears stabilized
        
        fprintf(' Chain appears stationary  but has "only"  %i samples  although the maximum of evaluations has been reached\n',(k-1-k_o));
        fprintf(' please re-run with longer MCMC_length  if you want a better sampling (check  histograms first)  \n');
        fprintf(fid_log,' Chain appears stationary  but has  "only"  ');fprintf(fid_log,num2str(k-1-k_o));fprintf(fid_log,' samples  although the maximum of evaluations has been reached\n');
        fprintf(fid_log,'please re-run with longer MCMC_length  if you want a better sampling (check  histograms first)  \n');
        
        
    end
    
end
warning on all;
 
Xf=X(1:k-1,:);
PXf=PX(1:k-1);

return
