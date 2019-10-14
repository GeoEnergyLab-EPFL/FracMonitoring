function [MPost,Ctilde_all,R_cor_all,W]=ProcessingMCMC(accept,PXf,Xf,k_o,stab,n_resampling,fpost)
%
%  perform diagnostic on MCMC results
%  
 
fprintf(' ----------------------- \n' );
fprintf(' PERFORM MCMC DIAGNOSTICS \n' );
fprintf(' Ellipse diffraction \n' );
fprintf(' ----------------------- \n' );

acceptrate=mean(accept);
fprintf('Acceptation rate: %f\n',acceptrate);
fprintf('Note that the acceptance rate should be between ideally 0.2 and 0.44\n');
fprintf(' - \n');
 

strings_m={'a','b','x_c','y_c','z_c','phi','theta','psi'};
col={'.r','.b','.m','.g','.y','.k','+r','+b'};

ns=n_resampling;
q=length(Xf(1,:)); % problem dimension

if  k_o==length(Xf(:,1))  
    k_o=k_o-1000;
end

tmp_ndx = find(accept, 1);

% resampling
if ( stab )
  k_s=[k_o:ns:length(PXf)];
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

disp(['Length of the chain after resampling:',num2str(length(k_s))]);

mLogP=PXf(k_s);
Xsp=Xf(k_s,:);

disp(['Fitting Gaussian multivariate pdfs mixture using BIC and EM_GM']);

[MPost,Ctilde_all,W]=fitGaussiansMixture(Xsp,fpost);

k=length(W);

for p=1:k

    CC=Ctilde_all(:,:,p);
     [Ctilde]=Posterior_Covariance_Matrix_ellipse(MPost(:,p));
    Ctilde_app(:,:,p)=Ctilde;
    %%% correlation matrix
    rho=Ctilde;
    R_cor=CC;
    for i=1:q
        for j=1:q
            rho(i,j)=Ctilde(i,j)/sqrt(Ctilde(i,i)*Ctilde(j,j));
            R_cor(i,j)=CC(i,j)/sqrt(CC(i,i)*CC(j,j));
        end
    end
    R_cor_all(:,:,p)=R_cor;
    rho_all(:,:,p)=rho;
end

%%%


%%% BY DEFAULT WE PLOT ONLY THE MOST PROBABLE MODE (the one corresponding
%%% with the highest weight.
CC=Ctilde_all(:,:,1);
R_cor=R_cor_all(:,:,1);
meanX=MPost(:,1);
Ctilde=Ctilde_app(:,:,1);
rho=rho_all(:,:,1);

handler_histm=figure('Name',' Histogram of model parameters from MCMC' );

% fit of Gaussian pdf to the mcmc results.

for i=1:q
    subplot(2,q/2,i);

%    x = linspace(min(Xsp(:,i)),max(Xsp(:,i)),100);
%    norm2 = normpdf(x,meanX(i),Ctilde(i,i)^0.5);
    hist(Xsp(:,i),60);hold on;
%    plot(x,norm2,'r-');
%    hold on; 
    
    if q==8
        title([' Hist of  ',strings_m{i}]); %, ' mu: ', num2str(mpost(i))
    else
        
        title([' Hist of m(',num2str(i),')']);
    end
    
end


theta=[0.:0.2:2*pi]';
xc=cos(theta);
yc=sin(theta);
myX=[xc yc];

 
handler_cor=figure('Name',' Correlation Plot between model Parameters with confidence ellipses ' ); %(from MCMC(green), from Quad. approx(magenta))
ng=(q-1);
kj=0;
% 
for i=1:q
    for j=i+1:q
        if j==i+1
            kj=kj+1+(i-1);
        else
            kj=kj+1;
        end
        subplot(ng,ng,kj);
        plot(Xsp(:,i),Xsp(:,j),'.');hold on;
        plot(meanX(i),meanX(j),'.r' ); hold on;
        
        %         % Draw confidence ellipse from Covariance matrix from chain
        Raux=[1  R_cor(i,j);
            R_cor(i,j)  1; ];
        
        aux=Raux*myX';
        aux(1,:)=aux(1,:)*2*sqrt(CC(i,i));   % 3 sigma
        aux(2,:)=aux(2,:)*2*sqrt(CC(j,j));
        
        plot(aux(1,:)+meanX(i) ,aux(2,:)+meanX(j),'g','LineWidth',1.5); hold on;
        
                 % Draw confidence ellipse from Covariance matrix from   Ctilde ::
        %
                 Raux=[1  rho(i,j);
                     rho(i,j)  1; ];
        %
                 aux=Raux*myX';
                 aux(1,:)=aux(1,:)*2.*sqrt(Ctilde(i,i));
                 aux(2,:)=aux(2,:)*2.*sqrt(Ctilde(j,j));
        %
                 plot(aux(1,:)+meanX(i) ,aux(2,:)+meanX(j),'m','LineWidth',1.5); hold on;
        %
        if q==8
            title([strings_m{i},'  vs  ',strings_m{j}]);
        else
            title(['m_',num2str(i),' vs m_',num2str(j)]);
        end
    end
end




%saveas(handler_histpost, 'HistogramOfLogPosterior.fig', 'fig');
%saveas(handler_histm, 'HistogramOfMpostFromMCMC.fig', 'fig');
%saveas(handler_cor, 'CorrelationPlotMpostConfidenceEllipses.fig','fig');
return

