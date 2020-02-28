%% This is an example of a way to load and look at acoustic data from a HF block experiment 
% GEL - EPFL - 2019
% Dong Liu -- 19/09/2019
% This script does not bring any difference in terms of the acoustic
% velocity change.
%
%% cleanup first, set global parameters
close all
clearvars
home

% data storage location 
% note here use 'gel-nas1' -  
% datastor = 'gel-nas1'; % 'local'; 'local' if dataset copied to local drive, 'gel-nas1', or 'enacdrives'
datastor = 'local';

% set the color, criteria, and change the size of the transducer to its
% original size

%% choose dataset and load acquisition times
switch datastor
    case 'gel-nas1'
        % data path on gel-nas1
        datapath = pathbyarchitecture('gel-nas1');
    case 'enacdrives'
        datapath = pathbyarchitecture('enac1files');
    case 'local'
        [~, username] = system('whoami');
        % datapath = ['/Users/' username(1:end-1) '/data/']; % here /Users is on OS X, on linux replace by home
        datapath = ['/Users/' username(1:end-1) '/Documents/Databackup/'];
end

% 2019 acquisitions

datayear = 19;

%test on gabbro
datamonth = 03;
dataday = 14;
starttime = '093621';

% test on marble
% datamonth = 07;
% dataday = 24;
% starttime = '095522';

% data folder name from experiment date
datafold = [num2str(datayear,'%02d') '-' num2str(datamonth,'%02d') '-' ...
    num2str(dataday,'%02d') '/'];

% extract header info from JSON file
fjson = [datapath datafold num2str(starttime) '.json'];
[jsonhdr,myTransducers,myPlattens,myBlock] = load_header(fjson);


% read active acoustic parameters from JSON header
ActiveAcoustic = jsonhdr.ActiveAcousticInfos;
% extract relevant parameters
np = ActiveAcoustic.NumberOfPoints;
ns = ActiveAcoustic.NumberOfSources;
nr = ActiveAcoustic.NumberOfReceivers;
Fs = ActiveAcoustic.SamplingFrequency_MHz_*1E6;
dt = 1/Fs;  % time step
t0 = 0;     % initial time
T = t0+dt*linspace(0,np-1,np)'; % time vector
Fn = 0.5*Fs;    % Nyquist frequency (Hz)

%%  read acoustic data from tge binary file 
fbin = [datapath datafold num2str(starttime) '.bin'];
seq_select=1:120;

%% load data of all pairs
[dataseq1, nq] = load_data(fbin,seq_select,1:32,1:32); 

%% load timing data
AcqTime = load_timing(fbin); % in date hours min sec format .... can be transformed in sec format

% time of the loaded sequence
AcSeqT=AcqTime(seq_select);

%% Set the velocity parameters to get the correct range for energy calculation
% for gabbro test
 Vp=6679;
 Vs=3668.5;
% for marble test
% Vp=5600;%6249.8;% to have better range for energy calculation
% Vs=3229.9;
endnoise=300; % burning period for EM noise from the source excitation
%%
transducerplot3D(myTransducers, myPlattens,[],[],1)
%% Quality Control of acoustic signal: SNR matrix on all pairs.
endnoise = 300; 
[normAll, h_SNR]  = signal_strength(squeeze(dataseq1(1,:,:,:)),endnoise)

%% MARB-005
i_pair=[17:22]';% for transducers north and south 5cm away from the platen center
% 17,sudden good signal quality after the fracture through the block
% 21,sudden disappearance of signal after the fracture through the block
% 20,22 very bad quality
% 18,19 shows good attenuation
i_pair=[31 32]'; % shear transducers north and south 2cm away from the platen center
% 32 shows good attenuation
% 31,sudden good signal quality after the fracture through the block 
i_pair=[23:28]';% for transducers west and east 2cm away from the platen center %24 is the shear
% 23,24(shear, but more obvious on the P-wave part) good attenuation of the signal
% 25,26 good attenuation but regain after the fracture going through
% 28 gain of energy after the fracture initiation
% 27 does not change that much
i_pair=[18 19 31 23 24 25 27]'
i_pair=[29 30]'% for shear transducers on the top and bottom plattens
%% GABB-001
i_pair=[17:22]'; % north and south 5cm away from the platen center % the change is within 5%
i_pair=[23 24 25 26 27 28]'; % east and west 2cm away from the platen center 
i_pair=[23 25 26]';
i_pair=24;
%24 is shear, 24,28 quality not good
% 24, 27 are in the center of the platen
% 17-22 no significant changes in amplitude within 5%
% we choose to plot Pair 23,25,26
%i_pair=[17 19 20 22 23 25 26 18 21 27]';% 24 28]';for gabbro
%i_pair=[29 30 31 32]';% shear transducers on the top and bottom plattens
i_pair=[23 25 26 27 30]';


%%
clear mywave_mode;
for i_d=1:length(i_pair)
    if myTransducers.wave_mode(i_pair(i_d))==1
        mywave_mode(i_d,:)='SS';
        %mywave_mode(i_d,:)='PP'; % we can force this to calculate the P-wave part
    else
        mywave_mode(i_d,:)='PP';
        %mywave_mode(i_d,:)='SS';% we can force in this way to calculate the S-wave part
    end
end

%SR_select = SourceReceiverPairs(myTransducers,myPlattens,[i_pair i_pair]);
SR_select = SourceReceiverPairs(myTransducers,myPlattens,[i_pair i_pair],mywave_mode);
fig_b = plotblockwithplattens(myBlock,myPlattens);
fig_handle = plotdirectrays(SR_select,fig_b);
%% loading the data

seqrange=seq_select;
[dataseq2, nq] = load_data(fbin,seqrange,i_pair',i_pair'); 
%[dataseq2, nq] = load_data(fbin,seqrange,i_pair,i_pair);

%% set the velocity range and plot the curve
tmid = floor((SR_select.distances/Vp)*Fs);
for i=1:SR_select.n_pairs
    if (SR_select.wave_type(i,1)=='S') % when the source is a shear transducer, we are more interested by its shear component
        tmid(i) = floor((SR_select.distances(i)/Vs)*Fs);
        disp(tmid(i));
    end   
end
windowsize=350;
tlow = tmid-windowsize;
tup = tmid+windowsize;   

% Plot the corresponding pair's waves
% set the plot range
tinmin=37*8000/160;%20
tinmax=45*8000/160;%160
% loop on each pair
for i=1:SR_select.n_pairs
    if i>1
        allwiggle = dataseq2(1:end,tinmin:tinmax,i,i);
    else
        allwiggle = dataseq2(1:end,tinmin:tinmax);
    end

dataSubswiggle = [zeros(1,size(allwiggle,2)); allwiggle(2:end,1:end)-allwiggle(1:end-1,1:end)];

%dataSubswiggle = dataSubswiggle - mean(dataSubswiggle,1);
%dataSubswiggle = dataSubswiggle - mean(dataSubswiggle,2);

% dataSubswiggle = allwiggle;
% dataSubswiggle = dataSubswiggle - mean(dataSubswiggle,2);

dataSubswiggle = [zeros(1,size(allwiggle,2)); allwiggle(2:end,1:end)-allwiggle(1,1:end)];


figure('DefaultAxesFontSize',24);
set(gcf,'NumberTitle','off');
set(gcf,'Name',['Transmission for S-R ' num2str(SR_select.SRmap(i,:))]);
%imagesc((1:size(allwiggle,1))',1e6*T(tinmin:tinmax)',(dataSubswiggle(1:end,1:end))',[min(min(dataSubswiggle)),max(max(dataSubswiggle))]);
wiggle(1e6*T(tinmin:tinmax)',(1:size(allwiggle,1))',(dataSubswiggle(1:end,1:end))','+');
colormap('gray')
axis ij
hold on
axis tight
hold on
yline(tlow(i)*160/8000,'-.r')
hold on
yline(tup(i)*160/8000,'-.r')
%line([pval pval], [1e6*T(endnoise) 1e6*T(end)],'Color', 'red','LineStyle', '--')
set(gcf,'Position',[100 100 500 700])
ylim([tinmin tinmax]*160/8000)
xlim([1 size(allwiggle,1)]) % this is the local sequence number
xlabel('Sequence number')
ylabel('Travel time (\mu s)')
end

%% pick the arrival if needed
splinedraw();
%% pick the arrival based on the spline SectionC5b
points= findobj(gca,'tag','trendpoints');
XData = get(points,'Xdata');
YData = get(points,'Ydata');
xmin=min(floor(XData+0.5));
xmax=max(floor(XData+0.5));
pickedTime=spline(XData,YData,xmin:1:xmax);
arrivalingg=[(xmin:1:xmax)' pickedTime'];

L_gg = plot(arrivalingg(1:end,1),arrivalingg(1:end,2),'k') % global plot results
set(L_gg,'ButtonDownFcn','Handle=gcbo');
disp('Do not close the figure and pick the arrival-time curve you want to save');

%% Section C7: get the data points from the chosen curve
arrival_time = Handle.YData;
sequence_local = Handle.XData;% what we pick here is the local sequence
sequence_picking = seq_select(sequence_local);

%% save the manual picking results into a txt file
close all
arrivalin = [sequence_picking' arrival_time'];
% please set your local saved path
% to be changed, once we need to save the data in the hf-gel-nas
localpath = 'SideTransmission/';
filename = ['S' num2str(SR_select.SRmap(1,1)) 'R' num2str(SR_select.SRmap(1,2)) SR_select.wave_type '.txt'];
seqnb = arrivalin(1:end,1); % here is the global sequence number

aux = [cellstr(datestr(AcqTime(seqnb))) num2cell(seqnb) num2cell(arrivalin(1:end,2))...
    num2cell((SR_select.SRmap(1,1)*ones(length(arrivalin),1))) num2cell(SR_select.SRmap(1,2)*ones(length(arrivalin),1))];
% we save the channel number here for source and receivers
% writecell function can be used directly in R2019a, instead of using a
% loop
fid = fopen([localpath filename],'wt');
for ii = 1:size(aux,1)% 
    fprintf(fid,'%s %.0f %f %.0f %.0f\n', aux{ii,:});
end
fclose(fid);


%% Calculate the energy

energy_all  = energy_calculation(dataseq2,endnoise,SR_select,1,[Vp Vs],Fs,350,1);

figure('DefaultAxesFontSize',24);
plot((1:length(energy_all))',energy_all)
% for gabbro
%legend('23','25','26')
%legend('23','24','25','26','27','28','23')
%legend('SR17','SR19','SR20','SR22','SR23','SR25','SR26','SR18','SR21','SR27')
%legend('29','30','31','32') % for shear transducers
legend('SR17','SR18','SR19','SR20','SR21','SR22')
% for marble
% legend('17','18','19','20','21','22')
% legend('31','32')
% legend('23','24','25','26','27','28')
% legend('18', '19', '31', '23', '24', '25', '27')
%legend('SR17','SR19','SR20','SR22','SR25','SR27','SR28')

%% save the attenuation signals into a txt file

%i_pair=[17 19 20 22 23 25 26 18 21 27]';% 10x1
energy_save=[i_pair'; energy_all];

fid = fopen('TransmissionSide.txt','wt');
for ii = 1:size(energy_save,1)% 
    fprintf(fid,'%g\t',energy_save(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);

%% plot transmitted signals on one plot
% similar to wiggle plot but different
i_pair=4;
tinmin=65*8000/160;
tinmax=75*8000/160;
i_switch=0;
figure
for i_seq=2:length(seq_select)
    hold on
    %plot(1e6*T(tinmin:tinmax)',dataseq2(i_seq,tinmin:tinmax,i_pair,i_pair)-0.005*i_switch)
    plot(1e6*T(tinmin:tinmax)',dataseq2(i_seq,tinmin:tinmax,i_pair,i_pair)-dataseq2(70,tinmin:tinmax,i_pair,i_pair)-0.005*i_switch)% removing the ref sequence
    i_switch=i_switch+1;
end
hold on
set(gcf,'Position',[100 100 500 700])
hold on
line([68 68], [-0.6 0],'Color', 'red','LineStyle', '--')%37.78 -7 0.05 for P-wave
%ylim([-12 2])
xlim([65 75]) %[35 45]
xlabel('Travel time (\mu s)')
ylabel('')
% S wave is more important

%% Try to pick the arrival one by one by inverse-plotting the signals

% set the local pair to plot
i_pair=1;
% set the time window to show
tinmin=20*8000/160;
tinmax=80*8000/160;
tendnoise=355*5;
pick_arrival=zeros(length(seq_select),1);
% set the switch 
i_switch=0;
figure
for i_seq=1:length(seq_select)
    hold on
    % plot the original data
    originalwave=dataseq2(i_seq,tinmin:tinmax,i_pair,i_pair);
    plot(1e6*T(tinmin:tinmax)',originalwave+i_switch*1)
    
    
    maxnoise=max(abs(dataseq2(i_seq,tinmin:tendnoise,i_pair,i_pair)));
    shiftwave=2*maxnoise-originalwave;
    diffwave=shiftwave-originalwave;
    id=find(diffwave<0);
    pick_arrival(i_seq)=id(1);
    hold on
    plot(1e6*T(tinmin:tinmax)',shiftwave+i_switch*1)
    i_switch=i_switch-1;
    
    hold on
    xlim([20 80]) %[35 45]
    xlabel('Travel time (\mu s)')
    ylabel('Amplitude (V)')
    % set the noise level
%     [~,noiselevel]=ginput(1);
%     % plot the shifted signal
%     plot(1e6*T(tinmin:tinmax)',-dataseq2(i_seq,tinmin:tinmax,i_pair,i_pair)+noiselevel)
%     % pick the arrival manually
%     [pick,~]=ginput(1);
%     pick_arrival(i_seq,1)=pick;
    hold on
    
end

% This method shows that the velocity increases with more damaged zones.
    %hold on
%set(gcf,'Position',[100 100 500 700])


% range for cross correlation can be 36.5 \mu s to 

% we can also plot the similar shifted plots for all



%% Velocity calculation for the selected pairs

%% plot transducer pairs to choose the period
% longitudinal transducers
figure
hold on
for i_pair = [1:23 25:28] %min([ns nr])
    plot(T*1E6,squeeze(dataseq1(15,:,i_pair,i_pair)))
end
xlabel('Time (\mus)')
ylabel('Amplitude (V)')
title('Longitudinal transducer pairs')

% shear transducers
figure
hold on
for i_pair = [24 29:32]
    plot(T(endnoise:8000)*1E6,squeeze(dataseq1(2,endnoise:8000,i_pair,i_pair)));
end
%legend('24','29','30','31','32')
xlabel('Time (\mus)')
ylabel('Amplitude (V)')
title('Shear transducer pairs')

%% calculate the velocity
% Separating P and S from class, create a matrix  S | R | distance | type | angle

%oppSR = transmitPairs;
%oppSR = AllPairsOppositePlattens(myTransducers,myPlattens); % change to not only opposite to not on the same platten !!! TO DO
oppSR=SR_select;
VmatrixIn.SRmap = oppSR.SRmap; % S/R
VmatrixIn.Ttype = string(oppSR.wave_type); % type P/S
VmatrixIn.Length = oppSR.distances; % length
VmatrixIn.N_pairs = oppSR.n_pairs; % n_pairs
VmatrixIn.Angle = oppSR.directions; % angle??
VmatrixIn.dt = dt; %dt, tipestep b/w measurements
 
%% load data for V only useful SR pair /// CREATE FOR ALL SEQ TBD!!!
for iSR=1:VmatrixIn.N_pairs % for all SRmap
    % VmatrixIn.dataset(:,iSR) = dataseq1(1,:,VmatrixIn.SRmap(iSR,1),VmatrixIn.SRmap(iSR,2)); % create a dataset in Vmatrix (dont comapre with D -- it is by receivers)
    VmatrixIn.dataset(:,iSR) = dataseq2(1,:,iSR,iSR);
end

figure
plot(T*1E6,VmatrixIn.dataset)
xlabel('Time (\mus)')
ylabel('Amplitude (V)')
%% 
% Ref
% dataset = Dd(:,1:5); % only opposite pairs and first seq ???
% reference signal read (without a sample), p-wave
VmatrixIn.dataRef_P = load_data([datapath '19-01-11/163814.bin'],1,1,1)'; % Ref signal:Swave = 162209, Pwave = 163814
VmatrixIn.dataRef_S = load_data([datapath '19-01-11/162209.bin'],1,1,1)'; % Ref signal:Swave = 162209, Pwave = 163814
% reference signal is only related to transducers, this is to remove the
% effect of the transducers
 
%% V velocity fuction
[VmatrixOut] = velocity_advanced(VmatrixIn, {35,45}, {60,90}, ["Full", "Part"]);%, "GI"]); % Vmatrix, {Pwave Tmin, Pwave Tmax}, {S wave Tmin, Swave Tmax}, ["Full", "Part", "GI", "AIC", "STA_LTA"], dt
%%
figure
plot(1:length(VmatrixOut.VpPart),VmatrixOut.VpPart(:))
figure
plot(1:length(VmatrixOut.VsPart),VmatrixOut.VsPart(:))

%% Calculate the velocity for several sequences
max_seq=90;
VpPart=zeros(max_seq,length(VmatrixOut.VpPart));
VsPart=zeros(max_seq,length(VmatrixOut.VsPart));
for i_seq=1:max_seq
    for iSR=1:VmatrixIn.N_pairs % for all SRmap
        %VmatrixIn.dataset(:,iSR) = dataseq1(i_seq,:,VmatrixIn.SRmap(iSR,1),VmatrixIn.SRmap(iSR,2)); % create a dataset in Vmatrix (dont comapre with D -- it is by receivers)
        VmatrixIn.dataset(:,iSR) = dataseq2(i_seq,:,iSR,iSR);
    end
    [VmatrixOut] = velocity_advanced(VmatrixIn, {35,45}, {65,75}, ["Full", "Part"]);
    VpPart(i_seq,:)=VmatrixOut.VpPart(:);
    %VsPart(i_seq,:)=VmatrixOut.VsPart([24 29:32]);
    % we calculate the Vs using the Vp signal by taking the S-arrival part,
    % the reference signal should be taken as the P-transducers signal,
    % since the reference signal serves as removing the transducer and
    % contact effect. 
    [VmatrixOut] = velocity_advanced(VmatrixIn, {35,45}, {65,75}, ["Full", "Part"]);
    VsPart(i_seq,:)=VmatrixOut.VsPart(:);
    close all;
end

% figure
% plot(1:16,VpPart(:,1:16)) % Vp=6249.8+-54.0
% figure
% plot(1:16,VsPart(:,1:16)) % Vs=3229.9+-176.2
%%
figure
plot(1:length(VpPart),VpPart(:,:))
% figure
% plot(1:length(VmatrixOut.VsPart),VmatrixOut.VsPart(:,:))

%%
figure
plot(1:length(VsPart),VsPart(:,:))
legend('SR23','SR25','SR26','SR27')

















