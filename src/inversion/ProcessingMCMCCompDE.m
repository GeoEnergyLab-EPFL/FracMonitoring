function [MPost,Ctilde_all,R_cor_all,W]=ProcessingMCMCCompDE(accept,PXf,Xf,k_o,stab,n_resampling,mpostde,sigpostde,fpost)
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
 


col={'.r','.b','.m','.g','.y','.k','+r','+b'};

ns=n_resampling;

%%%%%%%%%%%%%%%%%%%%29/05/2020%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q=length(Xf(1,:))-1; % problem dimension
clr3=[92 201 99]/255.;
clr4=[37 144 255]/255.;
clr0=[105 105 105]/255.;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global m_ind;

switch q
    case 8
        strings_m={'ln a','ln b','x_1','x_2','x_3','\psi','\theta','\phi'};
    case 7
        strings_m={'ln a','ln b','x_1','x_2','\psi','\theta','\phi'};
    case 6
        if(m_ind==2)
            strings_m={'ln r','x_1','x_2','x_3','\theta','\phi'};
        else
            strings_m={'ln a','ln b','x_1','x_2','x_3','\psi'};
        end
    case 5
        strings_m={'ln r','x_1','x_2','\theta','\phi'};
    case 4
        strings_m={'ln r','x_1','x_2','x_3'};
end


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


%%% the results from the DE algorithm
meanXD=mpostde;
[CtildeD]=Posterior_Covariance_Matrix_ellipse(mpostde);
rhoD=CtildeD;
for i=1:q
    for j=1:q
            rhoD(i,j)=CtildeD(i,j)/sqrt(CtildeD(i,i)*CtildeD(j,j));
    end
end



handler_histm=figure('Name',' Histogram of model parameters from MCMC' ,'DefaultAxesFontSize',18);
%set(gcf,'Position',[100 100 900 600]); % set the figure size;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 16, 8], 'PaperUnits', 'Inches', 'PaperSize', [16, 8])

% fit of Gaussian pdf to the mcmc results.

for i=1:q
    subplot(2,q/2,i); % this is only suitable for the case of 8 and 6.

    x = linspace(min(Xsp(:,i)),max(Xsp(:,i)),100);
    norm2 = normpdf(x,meanX(i),Ctilde(i,i)^0.5);
    [Nh,Xh]=hist(Xsp(:,i),60);
    Bh=bar(Xh,Nh,'facecolor',clr3);
    hold on;
    plot(x,norm2,'r-');
    hold on; 
    
    if q>=4 && q<=8
            title([' Hist of  ',strings_m{i}]); %, ' mu: ', num2str(mpost(i))    
    else
            title([' Hist of m(',num2str(i),')']);
    end
    
end


theta=[0.:0.2:2*pi]';
xc=cos(theta);
yc=sin(theta);
myX=[xc yc];

 
handler_cor=figure('Name',' Correlation Plot between model Parameters with confidence ellipses ' ,'DefaultAxesFontSize',14);
%set(gcf,'Position',[100 100 900 600]); % set the figure size;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 16, 8], 'PaperUnits', 'Inches', 'PaperSize', [16, 8])
%(from MCMC(green,black), from Quad. approx(magenta,blue))

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
        pointcolor=plot(Xsp(:,i),Xsp(:,j),'.');hold on;
        set(pointcolor,'color',clr3);
        plot(meanX(i),meanX(j),'.k' ); hold on;
        
        %         % Draw confidence ellipse from Covariance matrix from chain
        Raux=[1  R_cor(i,j);
            R_cor(i,j)  1; ];
        
        aux=Raux*myX';
        aux(1,:)=aux(1,:)*2*sqrt(CC(i,i));   % 3 sigma
        aux(2,:)=aux(2,:)*2*sqrt(CC(j,j));
        
        plot(aux(1,:)+meanX(i) ,aux(2,:)+meanX(j),'Color',clr0,'LineWidth',1.5); hold on;
        
                 % Draw confidence ellipse from Covariance matrix from   Ctilde from DE algorithm::
        %
                 Raux=[1  rhoD(i,j);
                     rhoD(i,j)  1; ];
        %
                 aux=Raux*myX';
                 aux(1,:)=aux(1,:)*2.*sqrt(CtildeD(i,i));
                 aux(2,:)=aux(2,:)*2.*sqrt(CtildeD(j,j));
                 
                 plot(meanXD(i),meanXD(j),'.b' ); hold on;
        %
                 plot(aux(1,:)+meanXD(i) ,aux(2,:)+meanXD(j),'Color',clr4,'LineWidth',1.5); hold on;
        %
        if q>=4 && q<=8
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

