function [MPost,Ctilde_all]=ProcessingMCMC(accept,PX,X,k_o,stab,mcmc_res_logfile)
%
% perform diagnostic on MCMC results
% Probably needs to re-direct the outputs to a kind of Log file
% Authour : Brice Lecampion
% date : May 27 2009
% date of last modification : July 9 2010
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
%
%  meanX  :: most probable value (mode) of the model param from MCMC
%          results
%  MPost(q,k) (q=6)  :: matrix containing the modes of the different
%  Gaussians pdf fitted on the MCMC sample
%
%  Ctilde_all(q,q,k) :: Posterior covariances around the mode  from quadratic
%          approximation (i.e Gaussian) for each fitted Gaussian


global G measured Cp prior_information


fprintf(' ----------------------- \n' );
fprintf(' PERFORM MCMC DIAGNOSTICS \n' );
fprintf(' ----------------------- \n' );

acceptrate=mean(accept);
fprintf('Acceptation rate: %f\n',acceptrate);
fprintf('Note that the acceptance rate should be between ideally 0.2 and 0.44\n');
fprintf(' - \n');

if (~isempty(mcmc_res_logfile))
    %fid_l=fopen(mcmc_res_logfile, 'wt');
    fid_l=fopen(mcmc_res_logfile, 'w');
    
    fprintf(fid_l,' ----------------------- \n' );
    fprintf(fid_l,'  MCMC DIAGNOSTICS \n' );
    fprintf(fid_l,' ----------------------- \n' );
    fprintf(fid_l,'Acceptation rate: %f\n',acceptrate);
    fprintf(fid_l,'Note that the acceptance rate should be between ideally 0.2 and 0.44\n');
    fprintf(fid_l,' - \n');
end

% resampling step
ns=3;
q=length(X(1,:)); % problem dimension
 
stab_length_min=500*q;  %3000 / 3 =   1000 samples  minimum left after re-sampling


global info_CP_save;

tmp_ndx = find(accept);
if ( isempty(tmp_ndx) )
  fprintf(1, '\n');
  fprintf(2, 'Warning: \n');
  fprintf(1, '	no stable result found after MCMC\n');
end
% X & PX only have accept results
% accept has 1/0 for all trials
%	numel(tmp_ndx) == numel(PX)

is_save = 1;


if (stab ~=  1 )
    %fprintf('Chain has not stabilized after %f accepted steps\n',length(PX));
    fprintf('Chain has not stabilized after %d accepted steps\n',length(PX));
    fprintf('The chain should be re-started (longer MCMCM_length)\n' );
    fprintf('No further processing will be done :: dubious results \n' );
    if (~isempty(mcmc_res_logfile))
        %fprintf(fid_l,'Chain has not stabilized after %f accepted steps\n',length(PX));
        fprintf(fid_l,'Chain has not stabilized after %d accepted steps\n',length(PX));
        fprintf(fid_l,'The chain should be re-started (longer MCMCM_length) \n' );
        fprintf(fid_l,'No further processing will be done :: dubious results \n' );
    end
    
    if length(PX)>1;
	fprintf(1, '\n');
	fprintf(2, 'Warning: \n');
	fprintf(1, '	only accepted MCMC used for statistics\n');
        MPost=mean(X(1:length(PX),:))';Ctilde_all=cov(X); %rho=[];CC=[];R_cor=[];
    else
        MPost=[];Ctilde_all=[];
    end

    if ( isfield(info_CP_save,'CP_save_mat') && info_CP_save.CP_save_mat >= 2 )
      fprintf(1, '\n');
      fprintf(2, 'Warning: \n');
      fprintf(1, '	info_CP_save.CP_save_mat = %d\n', info_CP_save.CP_save_mat);
      fprintf(1, 'generate MCMC histogram (no resampling) in any case\n');
      MCMC_hist(accept,PX,X,k_o,stab, is_save);
    end

    return;
    % the possibility of re-starting the chain IS NOT yet coded
    
elseif  (length(PX)-k_o) < stab_length_min
    
    %fprintf('Chain has stabilized after %f accepted steps\n',length(PX));
    fprintf('Chain has stabilized after %i accepted steps\n',k_o);
    fprintf(' ... and has run %i  steps (accepted) after that step \n' ,length(PX)-k_o);
    fprintf(' but it is not sufficiently sampled \n');
    fprintf(' The best is to increase MCMC_length and re-run the inversion \n' );
    fprintf('No further processing will be done at this stage. \n' );
    
    if (~isempty(mcmc_res_logfile))
        
        %fprintf(fid_l,'Chain has stabilized after %f accepted steps\n',length(PX));
        fprintf(fid_l,'Chain has stabilized after %i accepted steps\n',k_o);
        fprintf(fid_l,' ... and has run %i steps (accepted) after that step \n' ,length(PX)-k_o);
        fprintf(fid_l,'  but it is not sufficiently sampled \n');
        fprintf(fid_l,' The best is to increase MCMC_length and re-run the inversion \n' );
        fprintf(fid_l,'No further processing will be done at this stage.  \n' );
        
    end
    
    if length(PX)>1;
	fprintf(1, '\n');
	fprintf(2, 'Warning: \n');
	fprintf(1, '	only accepted MCMC used for statistics\n');
        MPost=mean(X(1:length(PX),:))';Ctilde_all=cov(X); %rho=[];CC=[];R_cor=[];
    else
        MPost=[];Ctilde_all=[];
    end

    if ( isfield(info_CP_save,'CP_save_mat') && info_CP_save.CP_save_mat >= 2 )
      fprintf(1, '\n');
      fprintf(2, 'Warning: \n');
      fprintf(1, '	info_CP_save.CP_save_mat = %d\n', info_CP_save.CP_save_mat);
      fprintf(1, 'generate MCMC histogram (no resampling) in any case\n');
      MCMC_hist(accept,PX,X,k_o,stab, is_save);
    end
    
    return;
    
else
    
    fprintf('Chain has stabilized after %i accepted steps\n',k_o);
    fprintf(' ... and has run %i  steps (accepted) after that step \n' ,length(PX)-k_o);
    fprintf('Processing of the results will be performed on the stationary part of the chain\n');
    fprintf(' Re-sampling the results every %i steps of the chain\n ',ns);
    if (~isempty(mcmc_res_logfile))
        fprintf(fid_l,'Chain has stabilized after %i accepted steps\n',k_o);
        fprintf(fid_l,' ... and has run %i steps (accepted) after that step \n' ,length(PX)-k_o);
        fprintf(fid_l,'Processing of the results will be performed on the stationary part of the chain\n');
        fprintf(fid_l,' Re-sampling the results every %i steps of the chain\n ',ns);
        
    end
    
end


% resampling
k_s=[k_o:ns:length(PX)];
 

mLogP=PX(k_s);
Xsp=X(k_s,:);

%%%% CHECK FOR MULTI_MODALITY using recursive EM
%%% FIT THE MCMC SAMPLING OF THE most probable parameters
%%% via a mixture of finite gaussian pdf via EM_GM algorithm
% we start with 5 gaussians, and decrease the number of gaussian until
% lowest BIC is found
% then we take the gaussian  with the highest weight as the MOST probable
% value

fprintf('---MULTI-MODES DETECTION ---\n');
fprintf('We now fit the MCMC sampling of the posterior pdf by a mixture of Gaussians \n');
fprintf('Starting with 6 Gaussians, and decreasing according to BIC criteria \n');
if (~isempty(mcmc_res_logfile))
    fprintf(fid_l,'---MULTI-MODES DETECTION ---\n');
    fprintf(fid_l,'We now fit the MCMC sampling of the posterior pdf by a mixture of Gaussians \n');
    fprintf(fid_l,'Starting with 6 Gaussians, and decreasing according to BIC criteria \n');
end

stfd=0; k=6; BIC_k_1=Inf;
while ~stfd && k>1
    
    [W_a,M_a,V_a,L,E] = EM_GM(Xsp,k,[],[],[],[])  ;
    BIC_k=-2*L+k*log(length(k_s));
    fprintf('BIC for %i Gaussians : %f \n',k,BIC_k);
    fprintf(fid_l,'BIC for %i Gaussians : %f \n',k,BIC_k);
    
    if  BIC_k<BIC_k_1 || norm(BIC_k/BIC_k_1-1) < 1e-3 % ENSURE SUFFICIENT INCREASE of BIC to stop with previous estimate
        
        k=k-1;
        BIC_k_1=BIC_k;
        W=W_a;M=M_a; V=V_a; 
        
    else

        k=k+1;
        stfd=1; % BIC re-increases, keep the previous results

    end
    
end

if k==1         % if exit with k=1, there is still 2 gaussians minimum in the mixture,
    k=k+1 ;
end
% sort the results in an ordered way, from the highest posterior value 
% this is where we particularize for function MinusLogPosterior etc.
for p=1:k,
        disp(M(:,p));
        m_log(p)=MinusLogPosterior(M(:,p));          
end

[m_log,IX]=sort(m_log,'ascend'); 
MPost=M(:,IX);
CC_all=V(:,:,IX);       % the associated fitted covariance matrix.
Worder=W(IX);


if k==2  % case of 2 pdfs fitted
    %%% checks that the 2 pdf are actually sufficiently distinct
    
    dis=max(abs((M(:,1)-M(:,2))./mean(M')'));  %L_inf norm of relative difference
    
    if dis < 0.5e-2  % the 2 modes are really close to one another (0.5%) : JUST MERGE
        
        MPost=M(:,1);
        CC_all=V(:,:,1);
        Worder=W*0+1;
        k=1; % only keep ONE Gaussian
    end
    
end

%%%% COMPUTE Covariance matrix, std deviation and correlation matrix for
%%%% each modes......
Ctilde_all=V;
for p=1:k

    [Ctilde]=Posterior_Covariance_Matrix(MPost(:,p),G, measured, Cp, prior_information);
    CC=CC_all(:,:,p);
    Ctilde_all(:,:,p)=Ctilde;
    %%% correlation matrix
    rho=Ctilde;R_cor=Ctilde;
    for i=1:q,
        for j=1:q,
            rho(i,j)=Ctilde(i,j)/sqrt(Ctilde(i,i)*Ctilde(j,j));
            R_cor(i,j)=CC(i,j)/sqrt(CC(i,i)*CC(j,j));
        end
    end
    R_cor_all(:,:,p)=R_cor;
    rho_all(:,:,p)=rho;
end

%%%% OUTPUTS to log file
fprintf(' %i  Gaussian(s) needed in this case \n',k);
if (~isempty(mcmc_res_logfile))
    fprintf(fid_l,' %f  Gaussians needed in this case\n',k);
end

%%% WE output all the most probables modes but in an ordered way:
for i=1:k
    fprintf('---------------------------------------------------------\n');
    fprintf('Properties of Gaussian number %i in the mixture \n',i);
    fprintf(' -Log Post. (at the mode) : %f \n', m_log(i));    
    fprintf(' Weight : %f \n', Worder(i));
    fprintf(' Most probable values of model parameters  = '); fprintf(num2str(MPost(:,i)'));fprintf(' \n');
    fprintf(' Post. Std Deviation (from MCMC)                 = ');fprintf(num2str( sqrt(diag(CC_all(:,:,i)))'));fprintf(' \n');
    fprintf(' Post. Std Deviation (from quad. app.)           = ');fprintf(num2str( sqrt(diag(Ctilde_all(:,:,i)))'));fprintf(' \n');
    fprintf(' Correlation matrix (MCMC): \n');
    cr=num2str(R_cor_all(:,:,i));
    for j=1:q
        fprintf(cr(j,:)); fprintf(' \n');
    end
    fprintf(' Correlation matrix (quad. app.): \n');
    cr=num2str(rho_all(:,:,i));
    for j=1:q
        fprintf(cr(j,:)); fprintf(' \n');
    end
    
    
    if (~isempty(mcmc_res_logfile))
        fprintf(fid_l,'---------------------------------------------------------\n');
        fprintf(fid_l,'Properties of Gaussian number %i in the mixture \n',i);
        fprintf(fid_l,' -Log Post. (at the mode) : %f \n', m_log(i));    
        fprintf(fid_l,' Weight : %f \n', Worder(i));
        fprintf(fid_l,' Most probable values of model parameters  = ');fprintf(fid_l,num2str(MPost(:,i)'));fprintf(fid_l,' \n');
        fprintf(fid_l,' Post. Std Deviation (from MCMC)                 = ');fprintf(fid_l,num2str( sqrt(diag(CC_all(:,:,i)))'));fprintf(fid_l,' \n');
        fprintf(fid_l,' Post. Std Deviation (from quad. app.)           = ');fprintf(fid_l,num2str( sqrt(diag(Ctilde_all(:,:,i)))'));fprintf(fid_l,' \n');
        fprintf(fid_l,' Correlation matrix (MCMC) : \n');
        cr=num2str(R_cor_all(:,:,i));
        for j=1:q
            fprintf(fid_l,cr(j,:)); fprintf(fid_l,' \n');
        end
        fprintf(fid_l,' Correlation matrix (quad. app.) : \n');
        cr=num2str(rho_all(:,:,i));
        for j=1:q
            fprintf(fid_l,cr(j,:)); fprintf(fid_l,' \n');
        end
        
    end
    
end
fprintf('---------------------------------\n');
fprintf('End of results from multi-modes detection \n' );

if (~isempty(mcmc_res_logfile))
    fprintf(fid_l,'---------------------------------\n');
    fprintf(fid_l,'End of results from multi-modes detection\n' );
end

%%% BY DEFAULT WE PLOT ONLY THE MOST PROBABLE MODE (the one corresponding
%%% with the highest weight.
CC=CC_all(:,:,1);
R_cor=R_cor_all(:,:,1);
meanX=MPost(:,1);
Ctilde=Ctilde_all(:,:,1);
rho=rho_all(:,:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOTS

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
    % plots the location of all the modes of the fitted gaussian
    for j=1:k
        plot(MPost(i,j),0,col{j}); hold on;
        %tmp_h{j} = plot(MPost(i,j),0,col{j}); hold on;
    end
    
    if q==6
        title([' Hist of  ',strings_m{i}]);
    else
        
        title([' Hist of m(',num2str(i),')']);
    end
end


%% to debug
%clear;
%openfig('HistogramOfMpostFromMCMC');
%k = 5;

% legend text
str_legend = '''hist(mpost)''';
for j = 1:k
  str_legend = sprintf('%s, ''# %d''', str_legend, j);
end
eval(sprintf('[h1, h2, h3, h4] = legend(%s);', str_legend));
%legend(h1, 'location', 'NorthEast');

% legend box
%Box_Pos = get(h1, 'Position')
%Box_Pos(1) = 1 - Box_Pos(3);
%Box_Pos(2) = 1 - Box_Pos(4);
%set(h1, 'Position', Box_Pos);

% inside legend: text/line/marker
if ( k > 3 )
  for i = 1:k
    %     hist
    itext = 1 + i;
    %     hist    hist
    iline = 1 + k + 1 + 2*i-1;
    imark = 1 + k + 1 + 2*i;

    %disp([i, itext, iline, imark]);

    % mark-X
    mark_XData = get(h2(imark), 'XData');
    set(h2(imark), 'XData', mark_XData/2);

    % text-X
    text_Pos = get(h2(itext), 'Position');
    text_Pos(1) = text_Pos(1) - mark_XData*3/4;
    set(h2(itext), 'Position', text_Pos);

    % line-X
    line_XData = get(h2(iline), 'XData');
    line_XData(2) = line_XData(2)/2;
    set(h2(iline), 'XData', line_XData);

    if ( i == 1 )
      YData_top = get(h2(imark), 'YData');
    elseif ( i == k )
      YData_bot = get(h2(imark), 'YData');
    end
  end

  YData_top = YData_top + (YData_top-YData_bot)/(k-1);
  YData_diff = (YData_top - YData_bot) / 3;
  for i = 1:3
    to_YData(i) = YData_top - YData_diff * i;
  end

  for i = 1:k
    % have the save Y
    % 	text (Position, 3D, XYZ)
    % 	line (YData, 2D, Y1Y2)
    % 	mark (YData, 1D)
    %i_to = i - 3;
    %imark_to = 1 + k + 1 + 2*i_to;
    %% mark
    %to_YData = get(h2(imark_to), 'YData');

    % id in h2
    %     hist
    itext = 1 + i;
    %     hist    hist
    iline = 1 + k + 1 + 2*i-1;
    imark = 1 + k + 1 + 2*i;

    if ( i <= 3 )
      i_toY = i;
      Xmov = 0;
    else
      i_toY = i - 3;
      Xmov = 0.5;
    end

    % mark-X
    mark_XData = get(h2(imark), 'XData');
    set(h2(imark), 'XData', mark_XData + Xmov);
    % mark-Y
    set(h2(imark), 'YData', to_YData(i_toY));

    % text
    text_Pos = get(h2(itext), 'Position');
    % text-X
    text_Pos(1) = text_Pos(1) + Xmov;
    % text-Y
    text_Pos(2) = to_YData(i_toY);
    set(h2(itext), 'Position', text_Pos);

    % line-X
    line_XData = get(h2(iline), 'XData');
    line_XData = line_XData + Xmov;
    set(h2(iline), 'XData', line_XData);
    % line-Y
    set(h2(iline), 'YData', [to_YData(i_toY), to_YData(i_toY)]);
  end
end



%%% CORRELATION Plots.....
handler_cor=figure('Name',' Correlation Plot between model Parameters with confidence ellipses ' ); %(from MCMC(green), from Quad. approx(magenta))
ng=(q-1);
k=0;

theta=[0.:0.2:2*pi]';
xc=cos(theta);
yc=sin(theta);
myX=[xc yc];

for i=1:q
    for j=i+1:q,
        if j==i+1
            k=k+1+(i-1);
        else
            k=k+1;
        end
        subplot(ng,ng,k);
        plot(Xsp(:,i),Xsp(:,j),'.');hold on;
        plot(meanX(i),meanX(j),'.r' ); hold on;
        
        %         % Draw confidence ellipse from Covariance matrix from chain
        Raux=[1  R_cor(i,j);
            R_cor(i,j)  1; ];
        
        aux=Raux*myX';
        aux(1,:)=aux(1,:)*2*sqrt(CC(i,i));   % 3 sigma
        aux(2,:)=aux(2,:)*2*sqrt(CC(j,j));
        
        plot(aux(1,:)+meanX(i) ,aux(2,:)+meanX(j),'g','LineWidth',1.5); hold on;
        
        %         % Draw confidence ellipse from Covariance matrix from   Ctilde ::
        %
        %         Raux=[1  rho(i,j);
        %             rho(i,j)  1; ];
        %
        %         aux=Raux*myX';
        %         aux(1,:)=aux(1,:)*1.*sqrt(Ctilde(i,i));
        %         aux(2,:)=aux(2,:)*1.*sqrt(Ctilde(j,j));
        %
        %         plot(aux(1,:)+meanX(i) ,aux(2,:)+meanX(j),'m','LineWidth',1.5); hold on;
        %
        if q==6
            title([strings_m{i},'  vs  ',strings_m{j}]);
        else
            title(['m_',num2str(i),' vs m_',num2str(j)]);
        end
    end
end

% %% SHOULD ask the user which ones he prefers Ctilde (quadratic
% %% approximation around the mode) or CC (direct MCMC results))

%%% WRITE ALL DIAGNOSTIC TO MCMC LOG FILE

if (~isempty(mcmc_res_logfile))
    
    cur_date=datestr(now);
    fprintf(fid_l, '-----MCMC Log file finished: \t: %s ----- \n\n', cur_date);
    fclose(fid_l);
    
end


saveas(handler_histpost, 'HistogramOfLogPosterior.fig', 'fig');
saveas(handler_histm, 'HistogramOfMpostFromMCMC.fig', 'fig');
saveas(handler_cor, 'CorrelationPlotMpostConfidenceEllipses.fig','fig');
return

