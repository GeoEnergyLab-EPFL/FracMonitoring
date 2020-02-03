function [mpost,sigpost]=MCMC_hist(accept,PX,X,k_o,stab,ns,varargin)

%
% function to generate MCMC histogram   
%  and fit independent gaussian (without correlation)  
% 
% INPUT PARAMETERS
%  accept :: rate of MCMC acception
%  X      :: MCMC vector outputs   a matrix number_of_MCMC_steps by number
%            of parameters
%  PX     :: a vector (length number_of_MCMC_steps) containing value
%            of -log Post
%  k_o    :: integer (define the end of the burn in period)
%  stab   :: integer flag (0, or 1) stabilized chain (1) or not (0)
%  ns     :: resampling step
%
% OUTPUTS
%	mpost :: mean of fitted individual normal pdf on each parameters
%   sigpost :: corresponding variance

q=length(X(1,:)); % problem dimension
global m_ind;


tmp_ndx = find(accept);

% resampling
if ( stab )
  k_s=[k_o:ns:length(PX)];
elseif ( ~isempty(tmp_ndx) )
  fprintf(1, '\n');
  fprintf(2, 'Warning: \n');
  fprintf(1, '	chain not stablized after MCMC\n');
  fprintf(1, '	only accepted results are used for plotting\n');
  k_s = 1:numel(PX);
else
  fprintf(1, '\n');
  fprintf(2, 'Warning: \n');
  fprintf(1, '	no stable result found after MCMC\n');
  fprintf(1, '	the last point is used for plotting\n');
  k_s = numel(PX);
end


mLogP=PX(k_s);
Xsp=X(k_s,:);

% save the Xsp into a txt file
fid = fopen('mLogP.txt','wt');
for ii = 1:size(Xsp,1)% 
    fprintf(fid,'%g\t',Xsp(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);



% -LOG POST histogram
handler_histpost=figure('Name','Histogram of - Log Posterior','DefaultAxesFontSize',24);
set(gcf,'Position',[100 100 900 600]); % set the figure size;

histfit(mLogP,60,'lognormal');  % it should look like a lognormal pdf
%title('Histogram of - Log Posterior');
%legend('It should look like a log Normal pdf' );


switch q
    case 8
        strings_m={'a','b','x_1','x_2','x_3','\psi','\theta','\phi'};
    case 7
        strings_m={'a','b','x_1','x_2','\psi','\theta','\phi'};
    case 6
        if(m_ind==2)
            strings_m={'r','x_1','x_2','x_3','\theta','\phi'};
        else
            strings_m={'a','b','x_1','x_2','x_3','\psi'};
        end
    case 5
        strings_m={'r','x_1','x_2','\theta','\phi'};
    case 4
        strings_m={'r','x_1','x_2','x_3'};
end

col={'.r','.b','.m','.g','.y','.k','+r','+b'};

handler_histm=figure('Name',' Histogram of model parameters from MCMC','DefaultAxesFontSize',24);
set(gcf,'Position',[100 100 900 600]) % set the figure size
% fit of Gaussian pdf to the mcmc results.
fitpostgaussian={};
mpost=zeros(q,1); 
sigpost=zeros(q,1); 
for i=1:q
    subplot(2,q/2,i);
    
    fitpostgaussian{i} = fitdist(Xsp(:,i),'Normal');
    mpost(i)= fitpostgaussian{i}.mu;
    sigpost(i)=fitpostgaussian{i}.sigma;

    histfit(Xsp(:,i),60);hold on;
    
    if q<=8 && q>=4
        title([' Hist of  ',strings_m{i}]); %, ' mu: ', num2str(mpost(i))
    else
        
        title([' Hist of m(',num2str(i),')']);
    end
    
end

 
return;
