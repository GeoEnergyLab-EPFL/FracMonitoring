% Inversion of real data
% using just fminsearch and approximate of posterior Covariance
%

%% cleanup first, set global parameters
close all
clearvars
home

% data storage location
datastor = 'local'; % 'local' if dataset copied to local drive, 'gel-nas1', or 'enacdrives'


%% choose dataset and load acquisition times
switch datastor
    case 'gel-nas1'
        % data path on gel-nas1
        datapath = pathbyarchitecture('gel-nas1');
    case 'enacdrives'
        datapath = pathbyarchitecture('enac1files');
    case 'local'
        [~, username] = system('whoami');
        %datapath = ['/home/' username(1:end-1) '/data/'];
        datapath = ['/Users/dongliu/AcousticHF/'];        
end

% 2019 acquisitions
datayear = 19;
% test on gabbro
datamonth = 03;
dataday = 14;
starttime = '093621';

% data folder name from experiment date
datafold = [num2str(datayear,'%02d') '-' num2str(datamonth,'%02d') '-' ...
    num2str(dataday,'%02d') '/'];

% extract header info from JSON file
fjson = [datapath datafold num2str(starttime) '.json'];
[jsonhdr,myTransducers,myPlattens,myBlock] = load_header(fjson);

%% Set experimental info

% diffraction data
% we read only the PP diffraction, one can change the wave_type by changing its value
wave_type = 'PP';
sidemarker = ['NV';'SV';'EV';'WV';'NN';'SS';'WW';'EE';'NE';'NW';'SE';'SW';'TT'];
fpath = [datapath datafold ];%'diffraction_picks/'];

%
% Set basic parameters for inverse problem
% gabbro = IsotropicSolid(3050,97.5867*1e9,0.3119); % vel in m/s before
% 28/11/2019
% relation between E, density and Poisson's ratio
gabbro = IsotropicSolid(3000,103.809*1e9,0.284314); % since 28/11/2019 

% set the input parameters
global d  SRPairs ray_type

sig_d = (0.5*1e-6);    % variance of measurement
% this should be the variance of the picked arrival time, can change from case to case
% we adopt here the average picked error 0.5\mu s

global Solid;
Solid = gabbro;

global prior

global m_ind % model indicator

%% M1_Ellipse model indicator as 1
%   a b    x y z, theta, phi psi
% guessed values vector a,b, center coordinate (XYZ),
% (alpha-local rotation; beta-tilt angle, most constrained; gamma-direction of the slope)
m_ind=1;
mp= [ log(0.05);log(0.05); .125;.125;.1285; 0.;0.;0. ];
sig_p = [-log(0.125);-log(0.125);0.02;0.02;0.006;pi/4;pi/60;pi/2]; % guessed variances
% relaxed tilt angle
%sig_p = [.125;.125;0.02;0.02;0.006;pi/4;pi/18;pi/2]; % guessed variances

%% M2_Radial model indicator as 2
% for the radial case r x y z alpha beta
m_ind=2;
mp= [ log(0.05);.125;0.125;.1285;0.;0.];
sig_p = [-log(0.125);0.02;0.02;0.006;pi/60;pi/2]; % guessed variances

%% M3_Ellipse with zero dip model indicator as 3
% for the elliptical case a b x y z self-rotation alpha
m_ind=3;
mp= [ log(0.05);log(0.05); .125;.125;.1285; 0.];
sig_p = [-log(0.125);-log(0.125);0.02;0.02;0.006;pi/4]; % guessed variances

%% M4_Radial with zero dip, model indicator as 4
% for the radial case with fixed z_c coordiante
m_ind=4;
mp= [ log(.05);.125;0.125;.1285];
sig_p = [-log(.125);0.02;0.02;0.006]; % guessed variances

%% M3bis_Ellipse with fixed z-coordinate for center model indicator as 5
% for the radial case r x y z alpha beta
m_ind=5;
global z_c
z_c=0.1285; % measured from the bottom
mp= [ 0.05;0.05;.125;.125;0.; 0.;0.];
sig_p = [.125;0.125;0.02;0.02;pi/4;pi/60;pi/2]; % guessed variances

%% M4bis_Radial with fixed z-coordinate for center model indicator as 6
% for the radial case with fixed z_c coordiante
m_ind=6;
global z_c
z_c=0.1285; % measured from the bottom
mp= [ .05;.125;.125;0.;0.];
sig_p = [.125;0.02;0.02;pi/60;pi/2]; % guessed variances

%% Set the prior
prior = GaussianPrior(mp,sig_p);% build the object

%% Set the calculate sequence range
seqrange=22:89;
%seqrange=28:30;
nseq=length(seqrange);

%% Inverse problem for M1-M4
% set the method to use
% value 0 indicates the quasi-Newton method
% value 1 indicates the Differential evolution method
switch m_ind
    case 1
% for M1
        inverse_method=[zeros(42,1);ones(nseq-42,1)];
% for M2
    case 2
        inverse_method=[zeros(49,1);ones(4,1);zeros(nseq-53,1)];
% for M3
    case 3
        inverse_method=[zeros(26,1);ones(nseq-26,1)];
    case 4
% for M4
        inverse_method=[zeros(57,1);1;zeros(nseq-58,1)];
    otherwise
        inverse_method=zeros(nseq,1);
end

inverse_method=ones(nseq,1);
inverse_method=zeros(nseq,1);

%seq40 seq79
mevol=zeros(nseq,length(mp)+1);% 22-114 %>=3traces:24-111 
sig_evol=zeros(nseq,length(mp)+1);
time={};
SRPairs_all={};
wrong_pick={};
i_w=0;
ray_type_all={};
prob_model=[];
Pair_number=[];% number of pairs usded for each sequence
likeli_model=[];

for i=1:nseq
    seqnb=seqrange(i);
    % Read the SRmap and arrival time
    [Pair_info, Pair_acqT]=load_diffraction(fpath, sidemarker, wave_type, seqnb);
    SRdiff = SourceReceiverPairs(myTransducers,myPlattens,Pair_info(:,1:2));
    Pair_number(i) = size(Pair_info,1);
    time{i}=Pair_acqT;
    SRPairs_all{i}= SRdiff;
    d = Pair_info(1:end,3)*1e-6; % arrival time data from diffraction should in s
    SRPairs = SRdiff; % SR pairs selected
    Cdinvdiag =(1/(1*sig_d^2))*ones(length(d),1); % data inverse of variance
    ray_type = ones(size(Pair_info,1),1);
    ray_type_all{i} = ray_type;
    
    % Check the S-R pairs
%     fig_b = plotblockwithplattens(myBlock,myPlattens)
%     fig_handle = plotdirectrays(SRdiff,fig_b);
    
    

    if inverse_method(i)==0
        %  direct unconstrained optimization (Quasi-Newton algorithm)
        %  set the MaxFunEvals
        options = optimset('MaxFunEvals',1e5);
        [z_opt,fval,exitflag,output] = fminsearch(@minusPosteriorPDF,[mp; log(sig_d)],options);
    else
        % differential evolution algorithm
        VTR = -1e6; % stop limit for the objective function, the value of the cost function below which the optimization would be stopped
        D = length(mp)+1;% number of parameters of the objective function
           
        % we don't give the intial input here, instead we set the minimum
        % and maximum values
        switch m_ind
            case 1
                XVmax = [log(.250);log(.250);.250;.250;.250;pi/4;pi/2;pi/2; log(20*sig_d)]'; % vector of upper bounds
                XVmin = [log(0.01);log(0.01);0.;0.;0.;-pi/4;-pi/2;-pi/2;log(0.01*sig_d)]';
            case 2
                XVmax = [log(.250);.250;.250;.250;pi/2;pi/2; log(20*sig_d)]';
                XVmin = [log(0.01);0.;0.;0.;-pi/2;-pi/2;log(0.01*sig_d)]';
            case 3
                XVmax = [log(.250);log(.250);.250;.250;.250;pi/2; log(20*sig_d)]'; % vector of upper bounds
                XVmin = [log(0.01);log(0.01);0.;0.;0.;-pi/2;log(0.01*sig_d)]';
            case 4
                XVmax = [log(.250);.250;.250;.250; log(20*sig_d)]';
                XVmin = [log(0.01);0.;0.;0.;log(0.01*sig_d)]';
            case 5
                % needs to write
                %XVmax = [log(.250);log(.250);.250;.250;pi;pi;pi/2; 10*sig_d]'; % vector of upper bounds
                %XVmin = [-5.;-5.;0.;0.;0.;0.;0.;-100]';
            case 6
                % needs to write
                %XVmax = [log(.250);.250;.250;pi;pi; 10*sig_d]';
                %XVmin = [-5.;0.;0.;0.;0.;-100]';
        end

        NP = 10*D; % number of population members, 10D is by defaut
        itermax = 2000; % maximum number of iterations (generations)
        F = 0.8; % DE-stepsize F from interval [0, 2]
        CR = 0.5; % crossover probability constant from interval [0, 1]
        strategy = 7;
        refresh = 10; 
        [z_opt,fval,nfeval] = devec3(@minusPosteriorPDF,VTR,D,XVmin,XVmax,[],NP,itermax,F,CR,strategy,refresh);
    end
    
    
    % quadratic approximation of the posterior variance
    [Cpost]=Posterior_Covariance_Matrix_ellipse(z_opt);
    sig_app_mpost=diag(Cpost).^0.5; % not sure about this.
    
    % calculate the Bayes factor
    cp=det(diag(sig_p.^2));
    % log version
    prob_model(i)=-fval-0.5*log(cp)+0.5*log(det(Cpost));
    % nonlog version
    %prob_model(i)=exp(-fval)/(cp.^0.5).*((det(Cpost)).^0.5);%/((2*pi).^((length(d)-1)/2));
   
    %%%%% check fit
    m = z_opt(1:length(z_opt)-1);
    noise = z_opt(length(z_opt));
    %SRPairs_multiseq{i_seq} = SRPairs;
    switch m_ind
        case 1
            ell = Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
        case 2
            ell = Radial(m(1),m(2:4),m(5),m(6));
        case 3
            ell = Ellipse(m(1),m(2),m(3:5),m(6),0,0);
        case 4
            ell = Radial(m(1),m(2:4),0,0);
        case 5
            ell = Ellipse(m(1),m(2),[m(3:4);z_c],m(5),m(6),m(7));
        case 6
            ell = Radial(m(1),[m(2:3);z_c],m(4),m(5));
        otherwise
            disp('Please check your input vector');
    end
    
    res = diffractionForward(Solid,SRPairs,ell,ray_type);% give one the shortest time needed for diffraction
    
    idx=find(abs(res(:,1)-d)*1e6 >10);
    if ~isempty(idx)
        i_w=i_w+1;
        wrong_pick{i_w,1}=seqnb;
        wrong_pick{i_w,2}=Pair_info(idx,1:2);
    end

%     % save the results in a json file 
     DiffRecord(i).seqnb=seqnb;
%     % acquisition time
     DiffRecord(i).acqT=time{i};
%     % SR_map used
     DiffRecord(i).mDE=m;
     % Model indicator
     DiffRecord(i).m_ind=m_ind;
     
    % calculate the likelihood
    mLikeli=minusLikelihoodEllipseDiff(z_opt);
    % better to use the log(likelihood) version
    likeli_model(i)=-mLikeli-0.5*length(d)*log(exp(noise)^2)-length(d)/2*log(2*pi);
 
    figure
    errorbar([1:length(d)]',d*1e6,ones(length(d),1)*exp(noise)*1e6,'b')
    hold on
    plot([1:length(d)]',res(:,1)*1e6,'*-r'); % the optimized arrival time
    xlabel('Source-Receiver Pair Number') % here the label is not clear it is the diffracted SR pick
    ylabel('Arrival Time (\mu s)')
    legend('real arrival time','calculated arrival time')
    hold on
    title(['Seq ' num2str(seqnb)])
    
    
    % write the info of the arrival and estimated error
%     if seqnb==50 % we choose the Seq.50
%         datasave=[[1:length(d)]',d*1e6,ones(length(d),1)*exp(noise)*1e6,res(:,1)*1e6];
%         fid = fopen('G01Seq50M2Matching.txt','wt');
%         for ii = 1:size(datasave,1)% 
%             fprintf(fid,'%g\t',datasave(ii,:));
%             fprintf(fid,'\n');
%         end
%         fclose(fid);
%     end
    
    
    mevol(i,:)=z_opt';
    sig_evol(i,:)=sig_app_mpost';
    
    % calculate the correlation matrix
    rho_matrix=ones(length(Cpost));
    for i_rho=1:length(Cpost)
        for j_rho=1:length(Cpost)
            rho_matrix(i_rho,j_rho)=sqrt(Cpost(i_rho,j_rho).^2/Cpost(i_rho,i_rho)/Cpost(j_rho,j_rho));
        end
    end
    % this is the correlation matrix containing log-parameters
end


%% Save the number of pairs used for each sequence
figure
plot(seqrange', Pair_number')

B=[seqrange' Pair_number'];
fileID = fopen(['GPair.txt'],'w');
fprintf(fileID,'%d %d\n',B');
fclose(fileID);

%% Save the correlation matrix
fileID = fopen(['G01Seq79Corr.txt'],'w'); % needs to adjust the name here
fprintf(fileID,'%d %d %d %d %d %d %d %d %d\n',rho_matrix);
fclose(fileID);

%% plot the correlation matrix
fid = fopen('G01Seq79Corr.txt'); % needs to adjust the name here
corr_matrix = textscan(fid,'%f%f%f%f%f%f%f%f%f','HeaderLines',0,'CollectOutput',1);
corr_matrix = corr_matrix{:};
fid = fclose(fid);

h=heatmap(abs(corr_matrix(1:8,1:8)))
colormap(viridis());

abscorr=abs(corr_matrix(1:8,1:8));

% fig_handle = figure('DefaultAxesFontSize',24);
% set(gcf,'Position',[100 100 600 600]);
fig_handle = figure('DefaultAxesFontSize',12);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3., 3.], 'PaperUnits', 'Inches', 'PaperSize', [3., 3.])
figure(fig_handle)
imagesc(1:size(abscorr,2),1:size(abscorr,1),abscorr)
axis square
caxis([0 1])
colormap(viridis());
colorbar
set(gca,'xtick',[]);
set(gca,'ytick',[]);

%% Adjust the data form
% save the real size, real estimated error and real Euler angles and
% modelling noise 
fracture_size=[];
off_set=[];
tilt_angles=[];
model_noise=[];
err_size_up=[];
err_size_low=[];
err_offset=[];
err_angles=[];
aspect_ratio=[];
err_ar_up=[];
err_ar_low=[];
    switch m_ind
        case 1 % elliptical
            fracture_size=exp(mevol(:,1:2)); % 2
            off_set=mevol(:,3:5); % 3
            tilt_angles=-1.*mevol(:,6:8); % 3
            model_noise=exp(mevol(:,9)); % 1
            
            err_size_up=exp(mevol(:,1:2)+sig_evol(:,1:2))-exp(mevol(:,1:2)); %2
            err_size_low=exp(mevol(:,1:2))-exp(mevol(:,1:2)-sig_evol(:,1:2)); %2
            err_offset=sig_evol(:,3:5); %3
            err_angles=sig_evol(:,6:8); %3
            
            aspect_ratio=exp(mevol(:,1)-mevol(:,2)); %1
            err_ar_up=exp(mevol(:,1)-mevol(:,2)+sqrt(sig_evol(:,1).^2+sig_evol(:,2).^2))...
                        -exp(mevol(:,1)-mevol(:,2)); %1
            err_ar_low=exp(mevol(:,1)-mevol(:,2))-exp(mevol(:,1)-mevol(:,2)-sqrt(sig_evol(:,1).^2+sig_evol(:,2).^2)); %1
                       
        case 2 % radial
            fracture_size=exp(mevol(:,1)); % 1
            off_set=mevol(:,2:4); % 3
            tilt_angles=-1.*mevol(:,5:6); % 2
            model_noise=exp(mevol(:,7)); % 1
            err_size_up=exp(mevol(:,1)+sig_evol(:,1))-exp(mevol(:,1)); %1
            err_size_low=exp(mevol(:,1))-exp(mevol(:,1)-sig_evol(:,1)); %1
            err_offset=sig_evol(:,2:4); %3
            err_angles=sig_evol(:,5:6); %2
        case 3 % elliptical
            fracture_size=exp(mevol(:,1:2)); % 2
            off_set=mevol(:,3:5); % 3
            tilt_angles=-1.*mevol(:,6); % 1
            model_noise=exp(mevol(:,7)); % 1
            err_size_up=exp(mevol(:,1:2)+sig_evol(:,1:2))-exp(mevol(:,1:2)); %2
            err_size_low=exp(mevol(:,1:2))-exp(mevol(:,1:2)-sig_evol(:,1:2)); %2
            err_offset=sig_evol(:,3:5); %3
            err_angles=sig_evol(:,6); % 1
            
            aspect_ratio=exp(mevol(:,1)-mevol(:,2)); %1
            err_ar_up=exp(mevol(:,1)-mevol(:,2)+sqrt(sig_evol(:,1).^2+sig_evol(:,2).^2))...
                        -exp(mevol(:,1)-mevol(:,2)); %1
            err_ar_low=exp(mevol(:,1)-mevol(:,2))-exp(mevol(:,1)-mevol(:,2)-sqrt(sig_evol(:,1).^2+sig_evol(:,2).^2)); %1
        case 4 % radial
            fracture_size=exp(mevol(:,1)); % 1
            off_set=mevol(:,2:4); % 3
            model_noise=exp(mevol(:,5)); % 1
            err_size_up=exp(mevol(:,1)+sig_evol(:,1))-exp(mevol(:,1)); %1
            err_size_low=exp(mevol(:,1))-exp(mevol(:,1)-sig_evol(:,1)); %1
            err_offset=sig_evol(:,2:4); %3
        case 5 % elliptical
            % needs to code
        case 6 % radial
            % needs to code
        otherwise
            disp('Please check your input vector');
    end

%% Save the results of inverse problems for M1
m_1=[fracture_size off_set tilt_angles model_noise aspect_ratio]; % 2+3+3+1+1
sig_1=[err_size_up err_size_low err_offset err_angles err_ar_up err_ar_low]; % 2+2+3+3+1+1
bayes_1=prob_model; % 1
wrongpick_1=wrong_pick; 
SRPairs_all_1=SRPairs_all;
ray_type_all_1=ray_type_all;
likelihood_1=likeli_model;
%% Save the results of inverse problems for M2
m_2=[fracture_size off_set tilt_angles model_noise aspect_ratio]; % 1+3+2+1
sig_2=[err_size_up err_size_low err_offset err_angles err_ar_up err_ar_low]; % 1+1+3+2
bayes_2=prob_model;
wrongpick_2=wrong_pick;
SRPairs_all_2=SRPairs_all;
ray_type_all_2=ray_type_all;
likelihood_2=likeli_model;
%% Save the results of inverse problems for M3
m_3=[fracture_size off_set tilt_angles model_noise aspect_ratio]; % 2+3+1+1+1
sig_3=[err_size_up err_size_low err_offset err_angles err_ar_up err_ar_low]; % 2+2+3+1+1+1
bayes_3=prob_model;
wrongpick_3=wrong_pick;
SRPairs_all_3=SRPairs_all;
ray_type_all_3=ray_type_all;
likelihood_3=likeli_model;
%% Save the results of inverse problems for M4
m_4=[fracture_size off_set tilt_angles model_noise aspect_ratio]; % 1+3+1
sig_4=[err_size_up err_size_low err_offset err_angles err_ar_up err_ar_low]; % 1+1+3
bayes_4=prob_model;
wrongpick_4=wrong_pick;
SRPairs_all_4=SRPairs_all;
ray_type_all_4=ray_type_all;
likelihood_4=likeli_model;


%% likelihood plot and modelling noise plot
figure
plot(likelihood_1)
hold on
plot(likelihood_2)
hold on
plot(likelihood_3)
hold on
plot(likelihood_4)
legend('M1','M2','M3','M4')

figure
plot(m_1(:,9))
hold on
plot(m_2(:,7))
hold on
plot(m_3(:,7))
hold on
plot(m_4(:,5))
legend('M1','M2','M3','M4')

%% Save into text files and post-process in mathematica
A1=[seqrange' m_1 sig_1 bayes_1' likelihood_1'];% 1+10+12+1+1
%%
A2=[seqrange' m_2 sig_2 bayes_2' likelihood_2'];% 1+7+7+1+1
%%
A3=[seqrange' m_3 sig_3 bayes_3' likelihood_3'];% 1+8+10+1+1
%%
A4=[seqrange' m_4 sig_4 bayes_4' likelihood_4'];% 1+5+5+1+1
%%
fileID = fopen(['G1.txt'],'w');
fprintf(fileID,'%d %f %f %f %f %f %f %f %f %e %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',A1');
fclose(fileID);
%%
fileID = fopen(['G2.txt'],'w');
fprintf(fileID,'%d %f %f %f %f %f %f %e %f %f %f %f %f %f %f %f %f\n',A2');
fclose(fileID);
%%
fileID = fopen(['G3.txt'],'w');
fprintf(fileID,'%d %f %f %f %f %f %f %e %f %f %f %f %f %f %f %f %f %f %f %f %f\n',A3');
fclose(fileID);
%%
fileID = fopen(['G4.txt'],'w');
fprintf(fileID,'%d %f %f %f %f %e %f %f %f %f %f %f %f\n',A4');
fclose(fileID);

%% Bayes factor plot
figure
semilogy(exp(bayes_2-bayes_1))
hold on
semilogy(exp(bayes_2-bayes_3))
hold on
semilogy(exp(bayes_2-bayes_4))
legend('B_{21}','B_{23}','B_{24}')

%% save the bayes factor
ABayes=[seqrange' (exp(bayes_2-bayes_1))' (exp(bayes_2-bayes_3))' (exp(bayes_2-bayes_4))' ];% 1+1+1+1

fileID = fopen(['GBayes.txt'],'w');
fprintf(fileID,'%d %e %e %e\n',ABayes');
fclose(fileID);

%% Check the quadratic approximation around the minimum for one certain m_post
% set the sequence number
seqnb=29;

% Read the SRmap and arrival time
[Pair_info, Pair_acqT]=load_diffraction(fpath, sidemarker, wave_type, seqnb);
SRdiff = SourceReceiverPairs(myTransducers,myPlattens,Pair_info(:,1:2));

d = Pair_info(1:end,3)*1e-6; % arrival time data from diffraction should in s
SRPairs = SRdiff; % SR pairs selected
Cdinvdiag =(1/(1*sig_d^2))*ones(length(d),1); % data inverse of variance
ray_type = ones(size(Pair_info,1),1);

[m_opt,fval,exitflag,output] = fminsearch(@minusPosteriorPDF,mp);

results=[];
for i=1:length(m_opt)
    for j=1:100
        m_new=m_opt;
        m_new(i)=m_opt(i)*(1+(j-50)*0.01); % from 0.5*m_opt to 1.5*m_opt
        results(j)=minusPosteriorPDF(m_new);
    end
    figure
    plot(m_opt(i)*(1+((1:100)'-50)./100),results')
    xlabel(['m' num2str(i)])
    ylabel('F*') % removed the constant terms
end

%% Save the output
% write into the json file 
txtoSave=jsonencode(DiffRecord);
fname='G01DForT.json';
fid=fopen(fname,'w');
fwrite(fid,txtoSave,'char');
fclose(fid);
% this is for transmission plot, saved sequence by sequence, different from
% that for FootPrint
% can be deduced from the FootPrint json file, needs to correct energy_plot
% function.

%% plot the evolution of the ellipse parameters
% set the estimated fracture initiation data
d2s=24*3600;
t_ini=3327/60; % should be adjusted with different tests
t_inj=datetime('14-Mar-2019 09:39:01');% 09:39:01
% set the time for the first sequence
t_first=(datenum(time{1})-datenum(t_inj))*d2s/60;

%% fracture size evolution

% what we get from the optimisation is acutually geometric mean and
% standard deviation of a=ln(A)

% we calculate the arithmetic mean and std of a .
% mevolnew=exp(mevol+0.5.*sig_evol.^2);
% sig_evolnew=exp(mevol+0.5.*sig_evol.^2).*sqrt(exp(sig_evol.^2)-1);

% we calculate the geometric mean and std of a.
mevolnew=exp(mevol);
%sig_evolnew=(exp(mevol+sig_evol)-exp(mevol-sig_evol))/2;
err_up = exp(mevol(1:nseq,1)+sig_evol(1:nseq,1))-mevolnew(1:nseq,1);
err_low = mevolnew(1:nseq,1)-exp(mevol(1:nseq,1)-sig_evol(1:nseq,1));


figure
% for one error bar
%errorbar([1:4:4*nseq],mevolnew(1:nseq,1),sig_evolnew(1:nseq,1))
% for two error bars: first lower error bar and then upper error bar
errorbar([1:4:4*nseq],mevolnew(1:nseq,1),err_low,err_up)
if (length(mp)==8)|| (length(mp)==7) || (m_ind==3)
    hold on
    errorbar([1:4:4*nseq],mevolnew(1:nseq,2),err_low,err_up,'r')
    legend('a','b')
    xlabel('Sequence number')
    ylabel('Elliptical radius (m)')
else
    legend('r')
    xlabel('Sequence number')
    ylabel('Radius (m)')
end
xticks(1:4*10:4*nseq)
xticklabels(seqrange(1:10:nseq))
%ylim([0 0.15])


% plot the aspect ratio
if (length(mp)==8)|| (length(mp)==7) || (m_ind==3)
    figure
    errorbar([1:4:4*nseq],aspect_ratio(1:nseq),err_ar_low(1:nseq),err_ar_up(1:nseq),'r')
    legend('M_1')
    xlabel('Sequence number')
    ylabel('Aspect Ratio')
end

%% fracture center coordiante evolution
figure
plot_idx=1; % radial case
if length(mp)==8 || length(mp)==7 || (m_ind==3)% ellipse case
    plot_idx=2;
end
errorbar([1:4:4*nseq],mevol(:,plot_idx+1),sig_evol(:,plot_idx+1))
hold on
errorbar([1:4:4*nseq],mevol(:,plot_idx+2),sig_evol(:,plot_idx+2),'r')
hold on
if (m_ind<=4)
    errorbar([1:4:4*nseq],mevol(:,plot_idx+3),sig_evol(:,plot_idx+3),'k')
end
legend('xc','yc','zc')
xticks(1:4*10:4*nseq)
xticklabels(seqrange(1:10:nseq))
xlabel('Sequence number')
ylabel('Centeral coordinate (m)')

%% Save the diffraction plot data from A1,A2,A3,A4
% we need to switch back the fracture size and tilting angles
model1_info=[log(m_1(:,1:2)) m_1(:,3:5) -m_1(:,6:8)];
%%
model2_info=[log(m_2(:,1)) m_2(:,2:4) -m_2(:,5:6)];
%%
model3_info=[log(m_3(:,1:2)) m_3(:,3:5) -m_3(:,6)];
%%
model4_info=[log(m_4(:,1)) m_4(:,2:4)];

%% save the results in a json file 
DiffRecord={};
DiffRecord.seqnb=seqrange';
% acquisition time
DiffRecord.acqT=time;

% Model indicator
DiffRecord.m_ind=2; % needs to change here
DiffRecord.mDE=model2_info; % needs to change here
DiffRecord.SRPairs=SRPairs_all_2; % needs to change here
DiffRecord.raytype=ray_type_all_2; % needs to change here

txtoSave=jsonencode(DiffRecord);
fname='G01ForFootPrintM2.json'; % needs to change here
fid=fopen(fname,'w');
fwrite(fid,txtoSave,'char');
fclose(fid);

%% read results from json code
% color styles
clr1=[68 1 84]/255.;
clr2=[253 231 37]/255.;
clr3=[92 201 99]/255.;
clr4=[37 144 255]/255.;
clr0=[105 105 105]/255.;

filename='/Users/dongliu/Documents/experimentDesignandResults/AcousticData/G01_14_03_2019/DiffractionInversion/SecondCheck/GNewInversion/G01ForFootPrintM1.json';
D1 = jsondecode(fileread(filename));
m_ind=D1.m_ind;
m_evol=D1.mDE;
SRPairs_all=D1.SRPairs;
ray_type_all=D1.raytype;
seqrange=D1.seqnb;

fig_handle = figure('DefaultAxesFontSize',14);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.75, 3.75], 'PaperUnits', 'Inches', 'PaperSize', [3.75, 3.75])
seqlist = [1:10:length(seqrange)];
seqrange(seqlist)
colorstyle=clr3;
% 2D fracture plot function
% input: sequence list, plot style
fig_handle=fractureFootprint(m_evol,seqlist,SRPairs_all,ray_type_all, myBlock,fig_handle,colorstyle,'-');
%%
hold on

filename='/Users/dongliu/Documents/experimentDesignandResults/AcousticData/G01_14_03_2019/DiffractionInversion/SecondCheck/GNewInversion/G01ForFootPrintM2.json';
D1 = jsondecode(fileread(filename));
m_ind=D1.m_ind;
m_evol=D1.mDE;
SRPairs_all=D1.SRPairs;
ray_type_all=D1.raytype;
seqrange=D1.seqnb;

colorstyle=clr0;
% 2D fracture plot function
% input: sequence list, plot style
fig_handle=fractureFootprint(m_evol,seqlist,SRPairs_all,ray_type_all, myBlock,fig_handle,colorstyle,'-');


%% Timing for the selected sequences
fbin = [datapath datafold num2str(starttime) '.bin'];
AcqTime = load_timing(fbin); % in date hours min sec format .... can be transformed in sec format
AcSeqT=AcqTime(seqrange);
tt=datenum(AcSeqT)-datenum(AcSeqT(1));
d2s=24*3600;
tt(seqlist)*d2s

%%
fig_b=plotblockwithplattens(myBlock,myPlattens);
fig_handle=plotdirectrays(SRdiff,fig_b);

%% output the data as a 3D video for one chosen model
% set the exception(errored) sequence
exception=[];
%fig2 = figure('units','normalized','outerposition',[0 0 1 1])
%fig2 = figure('DefaultAxesFontSize',24);
fig2 = figure('DefaultAxesFontSize',16);
%set(gcf,'Position',[100 100 600 600]); % set the figure size;
set(gcf,'Position',[100 100 500 500]); % set the figure size;
i_output = 0;
for i = 1:length(seqrange) %1:nseq
    m_i = m_evol(i,:);
    [~,id]=ismember(i,exception);
    if id==0
        fractureShape_plot(m_i,Solid,SRPairs_all(i),ray_type_all{i},...
            myBlock,myTransducers,myPlattens,fig2,'k.-');
        i_output=i_output+1;
        hold on
        title(['\fontsize{20}GABB-001 Seq ' num2str(seqrange(i)) ': ' datestr(AcqTime(i))])
        hold on
        text(-0.,-0.065,0,sprintf('� EPFL, Geo-Energy Lab (2019-2020) \nD. Liu, B. Lecampion, T. Blum'))
        F(i_output)=getframe(fig2);
            %pause(2)
        if i<nseq
            clf;
        end
    end

end

% create the video writer with 1 fps
  writerObj = VideoWriter('GABB-001-M1.mp4','MPEG-4');
  writerObj.FrameRate = 10;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);

%% analyse combined with the opening for one chosen model
% if considered as a radial fracture
% calculate the equivalent fracture radius
radius_eq=mevol(:,1);
if (length(m)==8)||(length(m)==7)
    radius_eq = sqrt(mevol(:,1).*mevol(:,2));
end
figure
plot(23+(1:82)',radius_eq)
xlabel('Sequence number')
ylabel('Fracture radius (m)')
% calculate the equivalent velocity
dura = zeros(length(time)-1,1);
for j=1:length(time)-1
    dura(j)=seconds(time{j+1}-time{j});
end
v_eq = (radius_eq(2:end)-radius_eq(1:end-1))./dura;
figure
plot(23+(1:82-1)',v_eq)
xlabel('Sequence number')
ylabel('Fracture front velocity (m/s)')
%% Adjust the velocity 
% get the indexes where the velocity is positive
idxpv=find(v_eq>0);
% one can smooth the v_calculation using moving average
%%
Ep=39.05.*1e9/(1-0.24^2);% we need check further for this value
Kp=2.84*1e6*sqrt(32/pi);
mup=0.6*12.;
ell_mk = Kp^6/Ep^4/mup^2./(v_eq(idxpv).^2) ;% this is very large
w_mk = Kp^4/mup./v_eq(idxpv)/Ep^3;
dimlsdist=radius_eq(idxpv+1)./ell_mk; % 10^-4-10^-5 the asymptote should follow the k-asymptote
%%
% plot the expected fracture opening and compare with the results for N-S
% and E-W
% we load the data from the SCCER-Soe Conference
widthPath = '/Users/dongliu/Documents/experimentDesignandResults/AcousticData/G01_14_03_2019/G01SCCERWidth/'; 
widthseqpath = [widthPath 'OpeningSequenceSCCERG01.txt'];
widthopeningpath = [widthPath 'OpeningSCCERG01.txt'];

% get the transmission coordinates
SRtrsn = SourceReceiverPairs(myTransducers,myPlattens,[1:16;1:16]');
trsn_x = SRtrsn.XS_XR(1:end,1); % x-coordinate
trsn_y = SRtrsn.XS_XR(1:end,2); % y-coordinate
trsn_z = 0.125*ones(16,1); % z-coordinate
% load the global sequence number
[seq_list] = importdata(widthseqpath,'\t');
% load the fracture opening
[width_profile] = importdata(widthopeningpath,'\t');

        
%% Opening profile plot
fig3 = figure('units','normalized','outerposition',[0 0 1 1]);
for j=1:82
    seq_i = 23+j; % the sequence number at which you want to plot the opening
    [~,idx] = ismember(seq_i,seq_list');

    if idx>0
        width_picked = (width_profile(idx,1:end))';
    else
        width_picked = zeros(16,1);
    end
    
    switch size(mevol,2)
        case 8
            c_x = mevol(j,3);
            c_y = mevol(j,4);
            c_z = mevol(j,5);
        case 7
            c_x = mevol(j,3);
            c_y = mevol(j,4);
            c_z = z_c;
        case 6
            c_x = mevol(j,2);
            c_y = mevol(j,3);
            c_z = mevol(j,4);
        case 5
            c_x = mevol(j,2);
            c_y = mevol(j,3);
            c_z = z_c;
    end
    w_eq = Kp/Ep*sqrt(radius_eq(j)).*sqrt(1-((1:100)/100).^2)*1e6;

    % N-S
    figure(fig3)
    widthplot1=fill([trsn_y(1);trsn_y(1:8); trsn_y(8)],[trsn_z(1);trsn_z(1:8)+width_picked(1:8);trsn_z(8)],'r');
    alpha(widthplot1,0.1)
    set(widthplot1,'EdgeColor','none')
    hold on
    title(['\fontsize{40}Seq ' num2str(j+24-1) ': ' datestr(time{j})])
    hold on
    plot((c_y-radius_eq(j)*(1:100)/100)' ,w_eq, 'r-')
    hold on
    plot((c_y+radius_eq(j)*(1:100)/100)' ,w_eq, 'r-' )
    
    % E-W
    hold on
    widthplot2 = fill([trsn_x(9); trsn_x(9:16); trsn_x(16)],[trsn_z(9);trsn_z(9:16)+width_picked(9:16);trsn_z(16)],'b');
    alpha(widthplot2,0.1)
    set(widthplot2,'EdgeColor','none')
    hold on
    plot((c_x-radius_eq(j)*(1:100)/100)' ,w_eq, 'b-')
    hold on
    plot((c_x+radius_eq(j)*(1:100)/100)' ,w_eq, 'b-' )
    xlim([0 myBlock.L_N])
    ylim([-10 60])
    legend('N-S Measurent','N-S toughness','estimation','E-W Measurement','E-W toughness','estimation')
%     F(j)=getframe(fig3);

    pause(0.1)
    if j<82
        clf;
    end
end


% % create the video writer with 1 fps
%   writerObj = VideoWriter('myVideoOpening.avi');
%   writerObj.FrameRate = 10;
%   % set the seconds per image
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:length(F)
%     % convert the image to a frame
%     frame = F(i) ;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);
%% Asymptote plot
fig4 = figure('units','normalized','outerposition',[0 0 1 1]);
loglog(10.^(-9:-2), (10.^(-9:-2)).^(1/2))
hold on
loglog(10.^(-4:0),2^(1/3)*3^(5/6)*(10.^(-4:0)).^(2/3))
xlabel('\xi_{mk}')
ylabel('\Omega_{mk}')

error_guess=0;
%error_guess=3e-6;

for j=1:82
    [~,id_j]=ismember(j,idxpv);
    if id_j>0
    seq_i = 23+j; % the sequence number at which you want to plot the opening
    [~,idx] = ismember(seq_i,seq_list');

    if idx>0
        width_picked = ((width_profile(idx,1:end))'.*1e-6+error_guess)./w_mk(id_j); % omega_mk
    else
        width_picked = zeros(16,1);
    end
    
    switch size(mevol,2)
        case 8
            c_x = mevol(j,3);
            c_y = mevol(j,4);
            c_z = mevol(j,5);
        case 7
            c_x = mevol(j,3);
            c_y = mevol(j,4);
            c_z = z_c;
        case 6
            c_x = mevol(j,2);
            c_y = mevol(j,3);
            c_z = mevol(j,4);
        case 5
            c_x = mevol(j,2);
            c_y = mevol(j,3);
            c_z = z_c;
    end
    ds = (radius_eq(j)- sqrt((trsn_y-c_y).^2+(trsn_x-c_x).^2))/ell_mk(id_j); % ksi_mk

    %plot
    figure(fig4)
    hold on
    title(['\fontsize{40}Seq ' num2str(j+24-1) ': ' datestr(time{j})])
    loglog(ds,width_picked,'.','MarkerSize',20);
%     loglog(ds(1),width_picked(1),'.','MarkerSize',20);
%     loglog(ds(2),width_picked(2),'*','MarkerSize',20);
%     loglog(ds(3),width_picked(3),'h','MarkerSize',20);
%     loglog(ds(4),width_picked(4),'+','MarkerSize',20);
%     loglog(ds(5),width_picked(5),'x','MarkerSize',20);
%     loglog(ds(6),width_picked(6),'s','MarkerSize',20);
%     loglog(ds(7),width_picked(7),'d','MarkerSize',20);
%     loglog(ds(8),width_picked(8),'p','MarkerSize',20);
%     loglog(ds(9),width_picked(9),'.','MarkerSize',20);
%     loglog(ds(10),width_picked(10),'*','MarkerSize',10);
%     loglog(ds(11),width_picked(11),'h-','MarkerSize',10);
%     loglog(ds(12),width_picked(12),'+-','MarkerSize',10);
%     loglog(ds(13),width_picked(13),'x-','MarkerSize',10);
%     loglog(ds(14),width_picked(14),'s-','MarkerSize',10);
%     loglog(ds(15),width_picked(15),'d-','MarkerSize',10);
%     loglog(ds(16),width_picked(16),'p-','MarkerSize',10);
    xlim([1e-9 1.])
    ylim([1e-5,10])

    hold on
    F(j)=getframe(fig4);
    end
end

% % create the video writer with 1 fps
%   writerObj = VideoWriter('G01Asymptote.avi');
%   writerObj.FrameRate = 1;
%   % set the seconds per image
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:length(F)
%     % convert the image to a frame
%     frame = F(i) ;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);

%% read the relection data
localpath = ['/Users/dongliu/AcousticHF/' datafold 'reflections/'];
filename = [num2str(datayear,'%02d') '-' num2str(datamonth,'%02d') '-' ...
    num2str(dataday,'%02d') '.txt'];

[arrival_info] = importdata([localpath filename],' ');
% this will import the data into textdata(date and time) and data (sequence number, arrival time and source number, receiver number)

source_channel = arrival_info.data(:,3); % this is channel+1 actually, for the sake of simplication we call the variable as channel
receiver_channel = arrival_info.data(:,4); 
picked_wavetype = char(arrival_info.textdata(:,1));
picked_seq1 = arrival_info.data(:,1);
picked_seq2 = arrival_info.data(:,5);

% plot the picked pair
picked_pairs = SourceReceiverPairs(myTransducers,myPlattens,[source_channel, receiver_channel],picked_wavetype);
fig_b = plotblockwithplattens(myBlock,myPlattens);
fig_handle = plotdirectrays(picked_pairs,fig_b);

%% plot of evolution together with reflection data
fig2 = figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:82
    m_i = mevol(i,:);
    fractureShape_plot(m_i,Solid,SRPairs_all{i},ray_type_all{i},...
        myBlock,myTransducers,myPlattens,fig2);
    hold on
    title(['\fontsize{40}Seq ' num2str(i+24-1) ': ' datestr(time{i})])
    hold on 
    idx_1=find(picked_seq1<=i+24-1);
    idx_2=find(picked_seq2<=i+24-1);

    % get the reflected points considering a flat fracture plane
    X_mid1=mean(picked_pairs.XS_XR(idx_1,[1 4]),2);
    Y_mid1=mean(picked_pairs.XS_XR(idx_1,[2 5]),2);
    X_mid2=mean(picked_pairs.XS_XR(idx_2,[1 4]),2);
    Y_mid2=mean(picked_pairs.XS_XR(idx_2,[2 5]),2);
    hold on
    plot3(X_mid1,Y_mid1,0.125*ones(length(X_mid1),1),'ro')
    hold on
    plot3(X_mid2,Y_mid2,0.125*ones(length(X_mid2),1),'bo')
    
    F(i)=getframe(fig2);
    %pause(2)
    if i<82
        clf;
    end
end
