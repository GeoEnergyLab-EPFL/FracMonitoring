%% This is an example of a way to load and look at acoustic data from a HF block experiment 
% GEL - EPFL - 2019
% Dong Liu -- 19/09/2019
% This script is better to be executed with a diffraction record json file.
% where .seqnb,  .mDE, and .acqT are least required.
%
%% cleanup first, set global parameters
close all
clearvars
home

% data storage location 
% note here use 'gel-nas1' -  
datastor = 'gel-nas1'; % 'local'; 'local' if dataset copied to local drive, 'gel-nas1', or 'enacdrives'

%% choose dataset and load acquisition times
switch datastor
    case 'gel-nas1'
        % data path on gel-nas1
        datapath = pathbyarchitecture('gel-nas1');
    case 'enacdrives'
        datapath = pathbyarchitecture('enac1files');
    case 'local'
        [~, username] = system('whoami');
        datapath = ['/Users/' username(1:end-1) '/data/']; % here /Users is on OS X, on linux replace by home
end

% 2019 acquisitions

datayear = 19;
%
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
[dataseq1, nq] = load_data(fbin,seq_select,1:32,1:32); 

%% load timing data
AcqTime = load_timing(fbin); % in date hours min sec format .... can be transformed in sec format

% time of the loaded sequence
AcSeqT=AcqTime(seq_select);

Vp=6679;
Vs=3534.1;
%% test the code for one-face plotting for the block
fig=plotside2Dwithplattens(myBlock,myPlattens,'B')
hold on
transducerplot2D(myTransducers,myPlattens,'B',fig,[],1)

%% Quality Control of acoustic signal: SNR matrix on all pairs.
endnoise = 300; % burning period for EM noise from the source excitation
[normAll, h_SNR]  = signal_strength(squeeze(dataseq1(1,:,:,:)),endnoise)

%% test the transmitted signals
[energy_all,Pairsmap]=energy_evolution(dataseq1,endnoise,myTransducers,myPlattens,[],[Vp,Vs],Fs,1/16*8000); % tranmistted signals

%% test the single platten singnals
[energy_all,Pairsmap]=energy_evolution(dataseq1,endnoise,myTransducers,myPlattens,15,'S',[],[Vp,Vs],Fs,1/16*8000); % single platten singles

%% Read the diffraction data
filename = 'DforT.json';
DiffractionRecord = jsondecode(fileread(filename));

%% set the energy structure and plot the results
Energy.seq=seq_select;
Energy.energy=energy_all;
Energy.SRmap=Pairsmap;
[F,fig_test]=energy_plot(Energy,myBlock,myPlattens,myTransducers,DiffractionRecord)

%% Save the results into a json file
% One can save the data after running all the sequences one wants to look at
% gloabl sequence  number

% write into the json file 
txtoSave=jsonencode(Energy);
fname='DTransmissionAllPartial.json';
fid=fopen(fname,'w');
fwrite(fid,txtoSave,'char');
fclose(fid);

%% Read the results from the saved json file
filename = 'myTforD.json';
Energy1 = jsondecode(fileread(filename));

fig_test=energy_plot(Energy1,myBlock,myPlattens,myTransducers,[])

%% Read the diffraction data
% plot the ellipse fracture and diffracted waves
fig2 = figure('units','normalized','outerposition',[0 0 1 1]);

%% write into a video file
% create the video writer with 1 fps
  writerObj = VideoWriter('TransmissionAllPartial.avi');
  writerObj.FrameRate = 4;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:100%length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


%% test the energy calculation code
oppositeSRs = SinglePlattenPairs(myTransducers,myPlattens,'T','B',2,'S');
energy_all  = energy_calculation(dataseq1,endnoise,oppositeSRs,[],[Vp,Vs],Fs,1/16*8000);
% one can comment the last four input arguments
figure
plot(energy_all);
