%% This is an example of picking the diffraction arrivals
% GEL - EPFL - 2019
% Dong Liu - 10/10/2019 

%% cleanup first, set global parameters
close all
clearvars
home

% data storage location 
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

% we use the gabbro experiment as an example
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

fbin = [datapath datafold num2str(starttime) '.bin'];

%% set the sequences you would like to read 
seqrange = 1:120;

%%  read acoustic data from the binary file 
[dataseq1, nq] = load_data(fbin,seqrange,1:32,1:32);
% continous sequences are much encouraged  

%% load timing data
AcqTime = load_timing(fbin); % in date hours min sec format .... can be transformed in sec format

%% loading the low rate data and set the estimated fracture initiation time
[lowrateData, lowrateTime, lowrateHeader] = load_lowrate(fbin);

% for better plot we need time in seconds (or minutes or hours) and 
% define the start of the injection as the initial zero time
% convert the datetime to seconds
d2s = 24*3600;
% the total injection rate
injectratevec = lowrateData(:,11)+lowrateData(:,14);
idtinject = find(injectratevec,1,'first');
% the datetime of the injection time
tinject = lowrateTime(idtinject); 
% time in min for the lowrate data
tlowrate = (datenum(lowrateTime)-datenum(tinject))*d2s/60; 
% time in seconds starting from the injection time for the selected sequences
tseq = (datenum(AcqTime(seqrange))-datenum(tinject))*d2s; 

% set the estimated initiation time as reference for the further
% diffraction picking
tinitiate = 3327; % should be adjusted with the time of fracture initiation
pval = tinitiate/60.;
% set the largest time of showing
tlowend = 200; % the largest shown time would be 200 min

% plot injection flows
figure
plot(tlowrate,lowrateData(:,11),tlowrate,lowrateData(:,14))
axis tight
xlim([tlowrate(1) tlowend])
ylim([0 1])
xlabel('Time (min)')
ylabel('Injection rate (mL/min)')

% plot pressure gauge and active acq times
figure
plot(tlowrate,lowrateData(:,2),tlowrate,lowrateData(:,3))
axis tight
xlim([tlowrate(1) tlowend])
ylim([0 40])
xlabel('Time (min)')
ylabel('Pressure (MPa)')
% add acq times
hold on
plot(tseq/60,ones(size(tseq))*10,'k.') % with unit as min
hold on
plot(tseq/60,ones(size(tseq))*10,'r.')
hold on
% add estimated initiation time
line([pval pval], [0 1.2*max(lowrateData(:,2))],'Color', 'red','LineStyle', '--')

%% Endnoise setting and Quality Control of acoustic signals
D = reshape(squeeze(dataseq1(1,:,:,:)),[],ns*nr); % flatten 3D array
Dd = squeeze(D(:,1:nr+1:end)); % extract 'diagonal'

figure
imagesc(1:nr,T*1E6,Dd)
xlabel('Receiver number')
ylabel('Time (\mus)')
title('All opposing pairs')
caxis([-1 1]*0.01)

endnoise = 250; % burning period for EM noise from the source excitation
[normAll, h_SNR]  = signal_strength(squeeze(dataseq1(1,:,:,:)),endnoise);


%% Diffraction starts from here.Choose the diffractions on one side N-S-E-W 
% During the whole picking process, please do not close the figures when 
% not asked to do so
% sidemarker = 'S'; % the input would be 'W','E','S'or'N'
% SRSidediff = SidePairsNeighborPlattens(myTransducers,myPlattens,sidemarker);
sidemarker = 'WT'; % WV->WT+WB
% SRSidediff = OnePlattenPairsAll(myTransducers,myPlattens,'N');->NN
SRSidediff = TwoPlattenPairsAll(myTransducers,myPlattens,'W','T')
fig_b = plotblockwithplattens(myBlock,myPlattens);
fig_handle = plotdirectrays(SRSidediff,fig_b);

%% To have a global look at all diffracted pairs from this side

[ival, iidex] = min(abs(tseq-tinitiate)); % find the nearest acquisition time
% draw the estimated initiation
if ival< tinitiate
    vmin = iidex;
    vmax = iidex+1;
else
    vmin = iidex-1;
    vmax = iidex;
end

pinterval = tseq(vmax)-tseq(vmin);
pval = vmin+(tinitiate-tseq(vmin))/pinterval;

% plot the diffraction curves for all the pairs
close all
[~] = diffraction_image(dataseq1, T, endnoise, SRSidediff, pval,0,[endnoise*dt np*dt]*10^6,11);
% 11 represents for the reference sequence, it can be also a range of
% sequences like [5 8 11]
% 1 represents that we plot the difference of the neighboring sequences, if
% it is 0, we plot the difference of all signals compared with the chosen
% reference sequences, for slowly propagating fracture, option 1 gives
% clearer diffracted curves
disp('Close the figures that you do not want to go further analysis');
%% Select the pairs with strong diffraction
figs = findobj('Type','figure');
fnb = zeros(length(figs),1);
for f = 1:length(figs)
    fnb1 = figs(f).Children(2).Title.String(4);
    fnb2 = figs(f).Children(2).Title.String(5);
    fnb3 = figs(f).Children(2).Title.String(6); % the three numbers.
    fnb(f,1) = str2num([fnb1 fnb2 fnb3]);
end
Selectmap = SRSidediff.SRmap(fnb,:);

% check the selected pairs
close all

SRs_select = SourceReceiverPairs(myTransducers,myPlattens,...
    SRSidediff.SRmap(fnb(:,1),:),SRSidediff.wave_type(fnb(:,1),:));
fig_b = plotblockwithplattens(myBlock,myPlattens);
fig_handle = plotdirectrays(SRs_select,fig_b);
disp(['You have ' num2str(length(fnb)) ' good pairs to look at on Side ' num2str(sidemarker)]);

%% to have a closer look with a given time range for this specific pair 
% and start to pick the time arrival
% for now it is better to deal with pair by pair since we can stop when
% there is a need % we chose the white wave valley point as the picked
% arrival
close all
pair_i = 2; % from 1 to length(Selectmap or fnb)
disp(['We are looking at S' num2str(SRSidediff.SRmap(fnb(pair_i,1),1)) 'R' num2str(SRSidediff.SRmap(fnb(pair_i,1),2))]);
SR_select = SourceReceiverPairs(myTransducers,myPlattens,SRSidediff.SRmap(fnb(pair_i,1),:),SRSidediff.wave_type(fnb(pair_i,1),:));

%initialize the arrival_time info
arrivalin = zeros(1,2);
arrivaling = arrivalin;
arrivalingg = arrivalin;

%% pick the arrival for each chosen sequence in the detailed plot where the chosen seq and only its neighboring sequences are shown in the plot
% plot the signals with one fixed reference signal is preferred
% one can leave the last argument empty to do this, if the last argument is
% not empty, one ask the function to plot the difference between the
% neighboring local sequences
close all
arrivalin = diffractions_picking(ActiveAcoustic,dataseq1, T, endnoise, SR_select, pval, 'local',[],10,[],[]);

%% pick the arrival based on the plot of the chosen sequences
close all
arrivaling = diffractions_picking(ActiveAcoustic,dataseq1, T, endnoise, SR_select, pval, 'semi-global',[],[],[],[]);

%% pick the arrival based on the plot of all the sequences loaded
close all
arrivalingg = diffractions_picking(ActiveAcoustic,dataseq1, T, endnoise, SR_select, pval, 'global',[],10,1,1);

%% check the arrival-picking results
close all
figshown = diffractions_picking(ActiveAcoustic,dataseq1, T, endnoise, SR_select, pval, ' ',[],11);
hold on 
L_g = plot(arrivaling(1:end,1),arrivaling(1:end,2),'r') % semi-global plot results
hold on 
L_gg = plot(arrivalingg(1:end,1),arrivalingg(1:end,2),'k') % global plot results
hold on 
L_l = plot(arrivalin(1:end,1),arrivalin(1:end,2),'b') % local plot results
legend('','estimated initiation','','semi-global','global','local')
set(L_g,'ButtonDownFcn','Handle=gcbo');
set(L_gg,'ButtonDownFcn','Handle=gcbo');
set(L_l,'ButtonDownFcn','Handle=gcbo');
disp('Do not close the figure and pick the arrival-time curve you want to save');

%% get the data points from the chosen curve
arrival_time = Handle.YData;
sequence_local = Handle.XData;% what we pick here is the local sequence
sequence_picking = seqrange(sequence_local);
%% write the arrival and the pair into a specific file
% we can write all into one file
% here for the sake of easy modification and re-picking for each pair, we
% save the pairs separately (at most we have ~100 pairs)
% the saving format is--the acquisition time, the local sequence number, 
% the source number and receiver number, and the arrival time

close all
arrivalin = [sequence_picking' arrival_time'];
% please set your local saved path
% to be changed, once we need to save the data in the hf-gel-nas
localpath = ['/Users/dongliu/AcousticHF/' datafold sidemarker '/'];
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

%% estimate the distance from the source to the notch and from the notch to the receiver
% we neglect the notch length here
s_xcor = SR_select.XS_XR(1:3);
r_xcor = SR_select.XS_XR(4:end);
notch_xcor = 0.5*myBlock.L_T*[1 1 1];% this can be improved by using the precise location
edge_xcor = 0.5*myBlock.L_T*[2 0 1];% this depends on the sidemarker 'N'
disdiff_first = vecnorm(s_xcor-notch_xcor)+vecnorm(r_xcor-notch_xcor);
mySolid.Vp = 6679;  
est_firstarrival = disdiff_first/(mySolid.Vp)*1e6; % in (\mus) we can calculate using the Vp for the first arrival

%% Read and check the picked arrival time from file for a specific pair
% data-loading takes a huge amount of the time, so the more economic way 
% is to load all the data in the beginning and check the chosen pairs
%
% set the path where you saved your arrival data
localpath = ['/Users/dongliu/AcousticHF/' datafold sidemarker '/'];
source_channel = 27; % got from the saved file name
receiver_channel = 9; % got from the saved file name
picked_wavetype = 'PP'; % got from the saved file name
filename=['S' num2str(source_channel) 'R' num2str(receiver_channel) picked_wavetype '.txt']; % should be better to create the side marker subfolder

% plot the picked pair
picked_pair = SourceReceiverPairs(myTransducers,myPlattens,[source_channel, receiver_channel],picked_wavetype);
fig_b = plotblockwithplattens(myBlock,myPlattens);
fig_handle = plotdirectrays(picked_pair,fig_b);

[arrival_info] = importdata([localpath filename],' ');
% this will import the data into textdata(date and time) and data (sequence number, arrival time and source number, receiver number)
picked_arrival = arrival_info.data(1:end,2); % in \mus

% try to re-check the picked arrival in the imagesc plot
picked_seq = string(arrival_info.textdata);
picked_seq = strcat(picked_seq(1:end,1)," ", picked_seq(1:end,2));
picked_date = datetime(picked_seq,'InputFormat','dd-MMM-yy HH:mm:ss');
% return sequence number 
[tf,picked_idx] = ismember(picked_date,AcqTime);% to find the absolute sequence number
if picked_idx == arrival_info.data(1:end,1)
    disp('Sequence member matching the time');
else
    disp('Sequence member does not match with the time');
end
% load the original data and plot 

endnoise = 250;
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % when we only want to check one pair or several pairs, we do not want%%%
% % % to load all the data                                                %%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% ref_seq=1; % set the reference signal, we take the first sequence here
% datacheck=load_data(fbin,[ref_seq picked_idx'],source_channel,receiver_channel);
% basesignal=datacheck(1,endnoise:end);
% allsignal=datacheck(1:end,endnoise:end);
% [dataSubs] = base_signal_substraction(allsignal,basesignal);
% 
% figure
% hold on
% title(['S' num2str(picked_pair.SRmap(1,1)) '-R' num2str(picked_pair.SRmap(1,2)) ' of type '  num2str(picked_pair.wave_type(1,:))  ' distance='  num2str(picked_pair.distances(1))],'fontsize',14);
% imagesc(1:(length(picked_idx)+1),1e6*T(endnoise:end),dataSubs',[min(min(dataSubs)),max(max(dataSubs))])
% colormap('gray')
% colorbar
% xticks((1:length(picked_idx)+1)')
% xticklabels([ref_seq picked_idx']')
% xlabel('Sequence number ()')
% ylabel('Travel time (\mu s)')
% hold on
% axis ij
% axis tight
% plot(2:length(picked_idx)+1,picked_arrival)

% % % % when we want to check all the pairs or more pairs let's say more than 5
% % % % pairs, it is better to load all the data in advance
% we need to find the local sequence in order to check
[~, local_idx] = ismember(picked_idx,seqrange);
close all
figcheck = diffractions_picking(ActiveAcoustic,dataseq1, T, endnoise, picked_pair, pval, ' ',[local_idx,picked_arrival],11,1,1);% We set Seq11 as the referenced sequence
title(['S' num2str(picked_pair.SRmap(1,1)) '-R' num2str(picked_pair.SRmap(1,2)) ' of type '  num2str(picked_pair.wave_type(1,:))  ' distance='  num2str(picked_pair.distances(1))],'fontsize',14);

%% Not satisfied? Correct arrivals for some sequences of this pair
% choose the sequences you wanna pick and rewrite the file with the corrected values
close all
figcorrect = diffractions_picking(ActiveAcoustic,dataseq1, T, endnoise, picked_pair, pval, ' ',[local_idx,picked_arrival],11,1);
title(['S' num2str(picked_pair.SRmap(1,1)) '-R' num2str(picked_pair.SRmap(1,2)) ' of type '  num2str(picked_pair.wave_type(1,:))  ' distance='  num2str(picked_pair.distances(1))],'fontsize',14);

disp('choose the sequences you wanna repick and finish selection with right click');
[seqnbc,] = getpts();% the last point is determined by the right click
seqnbc = floor(seqnbc+0.5); % the sequence number you wanna correct the arrival time
seqnbc = (sort(seqnbc))'; 
[~, correct_idx] = ismember(seqnbc,local_idx);
newpicked_arrival = picked_arrival;

for c=1:length(seqnbc)
    disp(['Pick the corrected arrival time for Seq' num2str(seqnbc(c))]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',['Correct arrival time for Seq' num2str(seqnbc(c))]);
    [~,correct_arrival] = ginput(1);
    newpicked_arrival(correct_idx(c),1) = correct_arrival;
    hold on
    plot(seqnbc(1:c),newpicked_arrival(correct_idx(1:c),1),'ro')
end

close all
[~] = diffractions_picking(ActiveAcoustic,dataseq1, T, endnoise, picked_pair, pval, ' ',[local_idx,picked_arrival],11);
title(['S' num2str(picked_pair.SRmap(1,1)) '-R' num2str(picked_pair.SRmap(1,2)) ' of type '  num2str(picked_pair.wave_type(1,:))  ' distance='  num2str(picked_pair.distances(1))],'fontsize',14);
hold on
plot(local_idx,newpicked_arrival)

%% rewrite all the info in the file with the same name
picked_arrival = newpicked_arrival;
% we always save the global sequence number
auxc =[cellstr(datestr(picked_date)) num2cell(picked_idx) num2cell(picked_arrival) num2cell((picked_pair.SRmap(1,1)*ones(length(picked_idx),1))) num2cell(picked_pair.SRmap(1,2)*ones(length(picked_idx),1))];
% we save the channel number here for source and receivers
% writecell function can be used directly in R2019a, instead of using a
% loop
fid = fopen([localpath filename],'wt');
for ii = 1:size(auxc,1)% 
    fprintf(fid,'%s %.0f %f %.0f %.0f\n', auxc{ii,:});
end
fclose(fid);

%% Bad results and I want to re-pick all the arrivals for this pair
% look at a specific pair and re-pick the arrival time from the beginning
% this is for the pair whose arrival time has been saved in .txt with a
% poor accuracy

% choose the arrival time gloablly 
SR_repick = SourceReceiverPairs(myTransducers,myPlattens,[source_channel, receiver_channel],picked_wavetype);
%initialize the arrival_time info

close all
arrival_repick = diffractions_picking(ActiveAcoustic,dataseq1, T, endnoise, SR_repick, pval, 'global',[local_idx, picked_arrival],11);% one can choose 'global', 'local' or 'semi-global'
seqnb = arrival_repick(1:end,1);
% we save always the global sequence number
seqnb = (seqrange(seqnb))';

auxrepick = [cellstr(datestr(AcqTime(seqnb))) num2cell(seqnb) num2cell(arrival_repick(1:end,2)) num2cell((SR_repick.SRmap(1,1)*ones(length(arrival_repick),1))) num2cell(SR_repick.SRmap(1,2)*ones(length(arrival_repick),1))];
%% rewrite the file with the same name
fid = fopen([localpath filename],'wt');
for ii = 1:size(auxrepick,1)% 
    fprintf(fid,'%s %.0f %f %.0f %.0f\n', auxrepick{ii,:});
end
fclose(fid);

