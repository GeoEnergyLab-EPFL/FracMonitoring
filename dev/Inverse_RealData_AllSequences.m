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
        datapath = ['/Users/bricelecampion/Documents/Work/Geomechanics/HydraulicFracturing/Experiments/local_data/'];
        
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

% load timing data
%fbin = [datapath datafold num2str(starttime) '.bin'];
%AcqTime = load_timing(fbin); % in date hours min sec format .... can be transformed in sec format

% %% Set the sequences that you would like to look at
% seq = [30 50 70 90];
% seq_n = length(seq);
% res_multiseq = zeros(seq_n,9); % sequence number + m-vector
% SRPairs_multiseq = cell(seq_n,1);
% % time of the loaded sequence
% %AcSeqT = AcqTime(seq);

%% Inversion on all sequences

% diffraction data
% we read only the PP diffraction, one can change the wave_type by changing its value
wave_type = 'PP';
sidemarker = ['N','S','E','W'];
fpath = [datapath datafold 'diffraction_picks/'];


%
% Set basic parameters for inverse problem
gabbro = IsotropicSolid(3050,97.5867*1e9,0.3119); % vel in m/s
% relation between E, density and Poisson's ratio

% set the input parameters
global d  SRPairs

sig_d = (0.5*1e-6);    % variance of measurement
% this should be the variance of the picked arrival time, can change from case to case
% we adopt here the average picked error 0.5\mu s

global Solid  Cdinvdiag;
Solid = gabbro;

global prior

%   a b    x y z, theta, phi psi
% guessed values vector a,b, center coordinate (XYZ), euler angle (alpha beta gamma)
mp= [ .05;.05; .125;.125;.125; 0.;0.;0. ];
sig_p = [.125;.125;0.02;0.02;0.005;pi/4;pi/40;pi/40]; % guessed variances


prior = GaussianPrior(mp,sig_p);% build the object

mevol=zeros(90-24,8);
sig_evol=zeros(90-24,8);
time={};
%seqnb =  90;%24
for i=1:(90-24)
    
    seqnb=23+i;
    % Read the SRmap and arrival time
    [Pair_info, Pair_acqT]=load_diffraction(fpath, sidemarker, wave_type, seqnb);
    SRdiff = SourceReceiverPairs(myTransducers,myPlattens,Pair_info(:,1:2));
    time{i}=Pair_acqT;
    d = Pair_info(1:end,3)*1e-6; % arrival time data from diffraction should in s
    SRPairs = SRdiff; % SR pairs selected
    Cdinvdiag =(1/(1*sig_d^2))*ones(length(d),1); % data inverse of variance
    
    % Check the S-R pairs
    % build the SRPair object
    % plot all the direct traces
%     fig_b = plotblockwithplattens(myBlock,myPlattens)
%     fig_handle = plotdirectrays(SRdiff,fig_b);
    
    %  direct unconstrained optimization (Nelder-Mead simplex)
    
    [m_opt,fval,exitflag,output] = fminsearch(@minusPosteriorPDF,mp);
    % quadratic approximation of the posterior variaance
    [Cpost]=Posterior_Covariance_Matrix_ellipse(m_opt);
    sig_app_mpost=diag(Cpost).^0.5;
    
    %%%%% check fit
    m = m_opt;
    %SRPairs_multiseq{i_seq} = SRPairs;
    ell = Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
    res = diffractionForward(Solid,SRPairs,ell);% give one the shortest time needed for diffraction
    
%     figure
%     errorbar([1:length(d)]',d*1e6,ones(length(d),1)*sig_d*1e6,'b')
%     hold on
%     %plot(d*1e6); hold on; % the real arrival time
%     %errorbar([1:length(d)]',res(:,1)*1e6,ones(length(d),1)*sig_d*1e6,'b')
%     plot([1:length(d)]',res(:,1)*1e6,'*-r'); % the optimized arrival time
%     xlabel('Source-Receiver Pair Number') % here the label is not clear it is the diffracted SR pick
%     ylabel('Arrival Time (\mu s)')
%     legend('real arrival time','calculated arrival time')
    
    mevol(i,:)=m_opt';
    sig_evol(i,:)=sig_app_mpost';
    
end

%%
figure
errorbar([1:4:4*66],mevol(:,1),sig_evol(:,1))
hold on
errorbar([1:4:4*66],mevol(:,2),sig_evol(:,2),'r')


figure
errorbar([1:4:4*66],mevol(:,3),sig_evol(:,3))
hold on
errorbar([1:4:4*66],mevol(:,4),sig_evol(:,4),'r')
hold on
errorbar([1:4:4*66],mevol(:,5),sig_evol(:,5),'k')
legend('xc','yc','zc')
