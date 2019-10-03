function MCMC_hist(accept,PX,X,k_o,stab, is_save)

%
% function to generate MCMC histogram for CP_save_mat >= 2 (full MCMC, no resampling), 
% in case it is not done in ProcessingMCMC (resampled with stride ns=3)
%
% INPUT PARAMETERS
%  accept :: rate of MCMC acception
%  X      :: MCMC vector outputs   a matrix number_of_MCMC_steps by number
%            of parameters
%  PX     :: a vector (length number_of_MCMC_steps) containing value
%            of -log Post
%  k_o    :: integer (define the end of the burn in period)
%  stab   :: integer flag (0, or 1) stabilized chain (1) or not (0)
%  mcmc_res_logfile
%         :: logfile if present
%
% OUTPUTS
%	N/A
%

if ( ~exist('is_save', 'var') || isempty(is_save) )
  is_save = 0;
end


% resampling step
ns=1;
q=length(X(1,:)); % problem dimension


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


% -LOG POST histogram
handler_histpost=figure('Name','Histogram of - Log Posterior');
hist(mLogP,60);  % it should look like a lognormal pdf
title('Histogram of - Log Posterior');
legend('It should look like a log Normal pdf' );


%%% Histogram plots of all the parameters

strings_m={'e_{xx}','e_{yy}','e_{xy}','f_{xx}','f_{yy}','f_{xy}'};
col={'.r','.b','.m','.g','.y','.k','+r'};

handler_histm=figure('Name',' Histogram of model parameters from MCMC' );
for i=1:q
    subplot(2,q/2,i);
    hist(Xsp(:,i),60);hold on;
    
    if q==6
        title([' Hist of  ',strings_m{i}]);
    else
        
        title([' Hist of m(',num2str(i),')']);
    end
end


fprintf(1, '\n');
fprintf(2, 'Warning: \n');
fprintf(1, '	skip Correlation Plot between model Parameters with confidence ellipses\n');

if ( is_save )
  fprintf(1, '	saving hist: log(pdf) & mpost\n');
  saveas(handler_histpost, 'CP3_HistogramOfLogPosterior.fig', 'fig');
  saveas(handler_histm, 'CP3_HistogramOfMpostFromMCMC.fig', 'fig');
else
  fprintf(1, '	skip saving figures\n');
end


%saveas(handler_histpost, 'HistogramOfLogPosterior.fig', 'fig');
%saveas(handler_histm, 'HistogramOfMpostFromMCMC.fig', 'fig');
%saveas(handler_cor, 'CorrelationPlotMpostConfidenceEllipses.fig','fig');

return;
