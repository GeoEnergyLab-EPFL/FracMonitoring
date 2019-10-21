% Inversion of real data
% using just fminsearch and approximate of posterior Covariance
%

%% cleanup first, set global parameters
close all
clearvars
home

% data storage location
datastor = 'gel-nas1'; % 'local' if dataset copied to local drive, 'gel-nas1', or 'enacdrives'


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
% mp= [ .05;.05; .125;.125;.125; 0.;0.;0. ];
% sig_p = [.125;.125;0.02;0.02;0.005;pi/4;pi/40;pi/40]; % guessed variances

% for the radial case r x y z alpha beta
mp= [ .05;.125;.125;.125; 0.;0.];
sig_p = [.125;0.02;0.02;0.005;pi/2;pi/40]; % guessed variances

prior = GaussianPrior(mp,sig_p);% build the object

mevol=zeros(111-24,length(m));% 22-114 %>=3traces:24-111 
sig_evol=zeros(111-24,length(m));
time={};
SRPairs_all={};
%seqnb =  90;%24

for i=1:(111-24)
    
    seqnb=23+i;
    % Read the SRmap and arrival time
    [Pair_info, Pair_acqT]=load_diffraction(fpath, sidemarker, wave_type, seqnb);
    SRdiff = SourceReceiverPairs(myTransducers,myPlattens,Pair_info(:,1:2));
    time{i}=Pair_acqT;
    SRPairs_all{i}= SRdiff;
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
    switch length(m)
        case 8
            ell = Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
        case 6
            ell = Radial(m(1),m(2:4),m(5),m(6));
        otherwise
            disp('Please check your input vector');
    end
    
    res = diffractionForward(Solid,SRPairs,ell);% give one the shortest time needed for diffraction

%     % save the results in a json file 
%     DiffRecord(i).seqnb=seqnb;
%     % acquisition time
%     DiffRecord(i).acqT=time{i};
%     % SR_map used
%     DiffRecord(i).mDE=m;
 
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
% write into the json file 
txtoSave=jsonencode(DiffRecord);
fname='DForT.json';
fid=fopen(fname,'w');
fwrite(fid,txtoSave,'char');
fclose(fid);


%% plot the evolution of the ellipse parameters
figure
errorbar([1:4:4*82],mevol(1:82,1),sig_evol(1:82,1))
if (length(m)==8)
    hold on
    errorbar([1:4:4*82],mevol(1:82,2),sig_evol(1:82,2),'r')
    legend('a','b')
else
    legend('r')
end

figure
plot_idx=1; % radial case
if length(m)==8 % ellipse case
    plot_idx=2;
end
errorbar([1:4:4*87],mevol(:,plot_idx+1),sig_evol(:,plot_idx+1))
hold on
errorbar([1:4:4*87],mevol(:,plot_idx+2),sig_evol(:,plot_idx+2),'r')
hold on
errorbar([1:4:4*87],mevol(:,plot_idx+3),sig_evol(:,plot_idx+3),'k')
legend('xc','yc','zc')

%% output the data as a video
fig2 = figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:82
    m_i = mevol(i,:);
    fractureShape_plot(m_i,Solid,SRPairs_all{i},...
        myBlock,myTransducers,myPlattens,fig2);
    hold on
    title(['\fontsize{40}Seq ' num2str(i+24-1) ': ' datestr(time{i})])
    F(i)=getframe(fig2);
    %pause(2)
    if i<82
        clf;
    end
end

% create the video writer with 1 fps
  writerObj = VideoWriter('myVideo.avi');
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

%% analyse combined with the opening
% if considered as a radial fracture
% calculate the equivalent fracture radius
radius_eq=mevol(:,1);
if length(m)==8
    radius_eq = sqrt(mevol(:,1).*mevol(:,2));
end
figure
plot(radius_eq)
xlabel('Sequence number')
ylabel('Fracture radius (m)')
% calculate the equivalent velocity
dura = zeros(length(time)-1,1);
for j=1:length(time)-1
    dura(j)=seconds(time{j+1}-time{j});
end
v_eq = (radius_eq(2:end)-radius_eq(1:end-1))./dura;
figure
plot(v_eq)
xlabel('Sequence number')
ylabel('Fracture front velocity (m/s)')
%%
Ep=3.*39.05.*1e9/(1-0.24^2);% we need check further for this value
Kp=2.84*1e6*sqrt(32/pi);
mup=0.6*12.;
ell_mk = Kp^6/Ep^4/mup^2./v_eq.^2 ;% this is very large
dimlsdist=radius_eq(2:end)./ell_mk; % 10^-4-10^-5 the asymptote should follow the k-asymptote
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
    
    c_x = mevol(j,3);
    c_y = mevol(j,4);
    c_z = mevol(j,5);
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
xlim([0.01 1.])
ylim([7,26])
loglog((1:80)/80, 20*((1:80)/80).^(1/2))
hold on
loglog((1:40)/40, 20*((1:40)/40).^(2/3))
xlabel('s/R')
ylabel('w (\mu m)')

for j=1:82
    seq_i = 23+j; % the sequence number at which you want to plot the opening
    [~,idx] = ismember(seq_i,seq_list');

    if idx>0
        width_picked = (width_profile(idx,1:end))';
    else
        width_picked = zeros(16,1);
    end
    
    c_x = mevol(j,3);
    c_y = mevol(j,4);
    c_z = mevol(j,5);
    ds = 1- sqrt((trsn_y-c_y).^2+(trsn_x-c_x).^2)./radius_eq(j);

    %plot
    figure(fig4)
    hold on
    title(['\fontsize{40}Seq ' num2str(j+24-1) ': ' datestr(time{j})])
    loglog(ds,width_picked,'.','MarkerSize',20);
    xlim([0.01 1.])
    ylim([7,26])

    hold on
    %F(j)=getframe(fig4);
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



