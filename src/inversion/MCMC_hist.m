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
% Arguments
% vararg1=mpost from DE algorithm
% vararg2=sigpost from DE algorithm

%%%%%%%%%%%%%%%%%% 29/05/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q=length(X(1,:))-1; % problem dimension, we remove the modelling noise as another dimension
% we can also add the modelling noise

clr3=[92 201 99]/255.;
clr0=[105 105 105]/255.;
clr4=[37 144 255]/255.;

narg = length(varargin);
if narg>=2
    if ~isempty(varargin)
        mpostde=varargin{1};
        sigpostde=varargin{2};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
set(gcf,'Units', 'Inches','Position',[0 0 8 6]); % set the figure size;

% Before we do histfit directly
% %hp=histfit(mLogP,60,'lognormal');  % it should look like a lognormal pdf
% 
% %hp=histfit(exp(-mLogP),60,'Normal');
% 
% hp=histogram(mLogP,60);
% %set(hp(1),'facecolor',clr3);
% set(hp,'Facecolor',clr3);
% %set(hp(2),'color',clr0);
% %title('Histogram of - Log Posterior');
% %legend('It should look like a log Normal pdf' );

% Now we first normalize the histogram and then do a fitting by ourselves
histogram(mLogP,'Normalization','pdf','FaceColor',clr3);
hold on
fitpostlog = fitdist(mLogP,'LogNormal');
mpostlog= fitpostlog.mu;
sigpostlog=fitpostlog.sigma;
xlog = linspace(min(mLogP),max(mLogP),100);
lognorm2 = lognpdf(xlog,mpostlog,sigpostlog);
plot(xlog,lognorm2,'Color',clr0,'Linewidth',2);





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

col={'.r','.b','.m','.g','.y','.k','+r','+b'};

handler_histm=figure('Name',' Histogram of model parameters from MCMC','DefaultAxesFontSize',18);
%set(gcf,'Position',[100 100 900 600]) % set the figure size
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 16, 8], 'PaperUnits', 'Inches', 'PaperSize', [16, 8])
% fit of Gaussian pdf to the mcmc results.
fitpostgaussian={};
mpost=zeros(q,1); 
sigpost=zeros(q,1); 
for i=1:q
    subplot(2,q/2,i);
    
    fitpostgaussian{i} = fitdist(Xsp(:,i),'Normal');
    mpost(i)= fitpostgaussian{i}.mu;
    sigpost(i)=fitpostgaussian{i}.sigma;

%     h=histfit(Xsp(:,i),60);
%     set(h(1),'facecolor',clr3);
%     set(h(2),'color',clr0);
%     hold on;
%     
%    if narg>=2
%     % get the coefficient in the plots
%     ymax=max(get(h(2),'YData'));
%     coef=ymax*sqrt(2*pi).*sigpost(i);
%     disp(coef);
    
    % plot the Gaussian distribution
%     x=get(h(2),'XData');
%     y=(1/sqrt(2*pi)/sigpostde(i)).*exp(-(x-mpostde(i)).^2/2/(sigpostde(i).^2));
%     %plot(x,y.*coef,'Color',clr4,'Linewidth',2);
% 
% % another way of plotting Gaussian distribution without adjusting the
% % coefficient
% %     x = linspace(min(Xsp(:,i)),max(Xsp(:,i)),100);
% %     norm2 = normpdf(x,mpostde(i),sigpostde(i));
% %     plot(x,norm2,'r-');
% %     hold on; 
%     hold on;
%     end
    
% arrange the histogram in a pdf way
    histogram(Xsp(:,i),'Normalization','pdf','FaceColor',clr3);
    hold on
    x = linspace(min(Xsp(:,i)),max(Xsp(:,i)),100);
    norm2 = normpdf(x,mpost(i),sigpost(i));
    plot(x,norm2,'Color',clr0,'Linewidth',2);
    
    y=(1/sqrt(2*pi)/sigpostde(i)).*exp(-(x-mpostde(i)).^2/2/(sigpostde(i).^2));
    plot(x,y,'Color',clr4,'Linewidth',2);
    hold on
    
    
    if q<=8 && q>=4
        title([strings_m{i}]); %, ' mu: ', num2str(mpost(i))
    else
        
        title([' Hist of m(',num2str(i),')']);
    end
    
end

 
return;
