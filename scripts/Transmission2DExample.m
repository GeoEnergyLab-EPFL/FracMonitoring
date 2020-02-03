%% This is an example of a way to load and look at acoustic data from a HF block experiment 
% GEL - EPFL - 2019
% Dong Liu -- 19/09/2019
% This script is better to be executed with a diffraction record json file.
% where .seqnb,  .mDE, and .acqT and .m_ind are least required.
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
seq_select=1:90;
[dataseq1, nq] = load_data(fbin,seq_select,1:32,1:32); 

%% load timing data
AcqTime = load_timing(fbin); % in date hours min sec format .... can be transformed in sec format

% time of the loaded sequence
AcSeqT=AcqTime(seq_select);

Vp=6679;
Vs=3668.5;
%% test the code for one-face plotting for the block
fig=plotside2Dwithplattens(myBlock,myPlattens,'B')
hold on
transducerplot2D(myTransducers,myPlattens,'B',fig,[],1)

%% Quality Control of acoustic signal: SNR matrix on all pairs.
endnoise = 300; % burning period for EM noise from the source excitation
[normAll, h_SNR]  = signal_strength(squeeze(dataseq1(1,:,:,:)),endnoise)

%% plot of the transmission signals for each pair
% plot all the pairs on the top and bottom
% detect all the transducers on the top and the bottom, regardless of the
% srouce or receiver 
for p=1:length(myPlattens)
    if myPlattens(p).face=='T'%'E'
        id_p=p;
    end
end
id_transducer=find(myTransducers.platten == myPlattens(id_p).id);
transmitnumber=myTransducers.channel(id_transducer)+1;
for i_d=1:length(id_transducer)
    if myTransducers.wave_mode(id_transducer(i_d))==1
        mywave_mode(i_d,:)='SS';
    else
        mywave_mode(i_d,:)='PP';
    end
end

% build the transducer Pairs and calculate the energy change
transmitPairs=SourceReceiverPairs(myTransducers,myPlattens,[transmitnumber transmitnumber],mywave_mode);
fig_b=plotblockwithplattens(myBlock,myPlattens);
fig_handle=plotdirectrays(transmitPairs,fig_b);

%% construct the diffraction pairs
oppositeSRs = SinglePlattenPairs(myTransducers,myPlattens,'T','B',8,'S');
fig_b=plotblockwithplattens(myBlock,myPlattens);
fig_handle=plotdirectrays(oppositeSRs,fig_b);

%% check the arrival of the diffraction/transmission plot
% give the Vp, Vs coresponding region
selectPairs=transmitPairs;
% selectPairs=oppositeSRs;
windowsize=350;
tmid = floor((selectPairs.distances/Vp)*Fs);
for i_t=1:selectPairs.n_pairs
    if (selectPairs.wave_type(i_t,:)=='SS') % when the source-receiver are shear transducers
        tmid(i_t) = floor((selectPairs.distances(i_t)/Vs)*Fs);
    end   
end
tlow = tmid-windowsize;
tup = tmid+windowsize;   

% check if the windowsize is good enough
for i_t=1:selectPairs.n_pairs
allwiggle = squeeze(dataseq1(1:end,endnoise:end,transmitnumber(i_t),transmitnumber(i_t)));
figure
set(gcf,'NumberTitle','off');
set(gcf,'Name',['Transmission for S-R ' num2str(i_t)]);
imagesc((1:size(allwiggle,1))',1e6*T(endnoise:end)',(allwiggle(1:end,1:end))',[min(min(allwiggle)),max(max(allwiggle))]);
colormap('gray')
axis ij
hold on
axis tight
hold on
yline(tlow(i_t)*160/8000,'-.r')
yline(tup(i_t)*160/8000,'-.r')
%line([pval pval], [1e6*T(endnoise) 1e6*T(end)],'Color', 'red','LineStyle', '--')
set(gcf,'Position',[100 100 500 700])
ylim([endnoise 8000]*160/8000)
xlim([1 size(allwiggle,1)]) % this is the local sequence number
xlabel('Sequence number')
ylabel('Travel time (\mu s)')
end

%% test the transmitted signals
[energy_all,Pairsmap]=energy_evolution(dataseq1,endnoise,myTransducers,myPlattens,[],[],16,[Vp,Vs],Fs,350); % tranmistted signals

%% test the single platten singnals
[energy_all,Pairsmap]=energy_evolution(dataseq1,endnoise,myTransducers,myPlattens,8,'S',11,[Vp,Vs],Fs,350); % single platten singles

%% plot the energy evolution with time
seq_select=1:90;
figure
for i=1:size(energy_all,2)
    energy_ratio=energy_all(:,i);
    hold on
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',['Transmission for S-R ' num2str(Pairsmap(i,:))]);
    if (i==19) || (i==20) || (i==9) || (i==10)
        plot(seq_select,energy_ratio')
    else
        %plot(seq_select,energy_ratio')
    end
    xlabel('Sequence number')
    ylabel('Energy ratio')
end
hold on 
legend('9-SR30', '10-SR32','19-SR29','20-SR31')
%% Read the diffraction data
filename = 'G01DforT.json';
DiffractionRecord = jsondecode(fileread(filename));

%% set the energy structure and plot the results
Energy.seq=seq_select;
Energy.energy=energy_all;
Energy.SRmap=Pairsmap;
[F,fig_test]=energy_plot(Energy,myBlock,myPlattens,myTransducers,DiffractionRecord,[],'b',[])

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

%% combined with the fracture opening to replace the energy
% should set the criteria differently
% we load the data from the SCCER-Soe Conference
% widthPath = '/Users/dongliu/Documents/experimentDesignandResults/AcousticData/G01_14_03_2019/G01SCCERWidth/'; 
% widthseqpath = [widthPath 'OpeningSequenceSCCERG01.txt'];
% widthopeningpath = [widthPath 'OpeningSCCERG01.txt'];
% new opening results
widthPath = '/Users/dongliu/Documents/experimentDesignandResults/AcousticData/G01_14_03_2019/GabbroWidth/Newresults/'; 
widthseqpath = [widthPath 'OpeningSequenceFullG01.txt'];
widthopeningpath = [widthPath 'OpeningFullG01.txt'];

% load the global sequence number
[seq_list] = importdata(widthseqpath,'\t');
% load the fracture opening
[width_profile] = importdata(widthopeningpath,'\t');

Energy.seq=seq_list;
Energy.energy=width_profile;
Energy.SRmap=[(1:16)' (1:16)'];
fig_test = figure('DefaultAxesFontSize',24);
set(gcf,'Position',[100 100 600 600]); % set the figure size;
[F,fig_test]=energy_plot(Energy,myBlock,myPlattens,myTransducers,DiffractionRecord,fig_test,'r',6)


%% plot the energy of one sequence % Sequence plot choosing from Seq.27, 28, 29, 30, 31, set criteria as 0.8
% we can set the plotted sequences % 27 30 31
%clf(fig_test,'reset');
sel_seq=75;%[22 27 31 39 65 ];% 85
Energy.seq=sel_seq;
% for transmission P-wave signals
Energy.energy=width_profile(sel_seq,:);
Energy.SRmap=[(1:16)' (1:16)'];
fig_test = figure('DefaultAxesFontSize',24);
set(gcf,'Position',[100 100 600 600]); % set the figure size;
[F,fig_test]=energy_plot(Energy,myBlock,myPlattens,myTransducers,DiffractionRecord,fig_test,'k',6)

%% plot only the 2D diffraction curve
sel_seq=58;
Energy.seq=sel_seq;
% for transmission P-wave signals
Energy.energy=ones(length(sel_seq),16);
Energy.SRmap=[(1:16)' (1:16)'];
fig_test = figure('DefaultAxesFontSize',24);
set(gcf,'Position',[100 100 600 600]); % set the figure size;
[F,fig_test]=energy_plot(Energy,myBlock,myPlattens,myTransducers,DiffractionRecord,fig_test,'k',0.8)


