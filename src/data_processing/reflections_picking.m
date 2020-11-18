function [arrival_info]=reflections_picking(activeInfo,dataseq1, T, endnoise, mySRPair, seqInitiate, plotoption, varargin)
% Dong Liu -- 31/10/2019
% the detailed plot for picking the arrival time for a certain pair
% in order to have a closer look for reflection data
% This function is using the function of diffraction_images.m, wiggle.m,
% base_signal_substraction.m
%
% activeInfo: active acoustic settings from the header, corresponding to
% the global variable ActiveAcoustic
% dataseq1: the loaded acoustic data
% mySRPair: chosen transducer pairs object
% seqInitiate: the guessed number of the sequence corresponding to the
% fracture initiation, obtained from the presssure responses
% diffindicator: bool value, if the value==1 plot the difference between the
% neighboring sequences
% plotoption: one only has two choices for plotoptions:
% global and semi-global
% optional argument1: give the (nx2) matrix (previously saved data for this pair), in
% case of the first-time picking, this matrix value is all 0--> zeros(1,2)
% optional argument2: give the reference sequence number 
% or the sequence ranges where the the averaged reference
% sequence is calculated from
% optional argument3: indicate whether to plot the wiggle plot or not for
% the plotting part
% optional argument4: indicate whether to plot the neighboring difference
% it is deactived by default
%
%% determine whether to plot the differences between neighboring sequences
narg = length(varargin);
neborplot=0; % the plot option is deactived by default.
if narg>=4 && ~isempty(varargin{4})
    disp('Plot only the difference between neighboring sequences')
    neborplot=1;
end

%% get the active acoustic info
np = activeInfo.NumberOfPoints;
Fs = activeInfo.SamplingFrequency_MHz_*1E6;
dt = 1/Fs;  % time step
nseq = size(dataseq1,1);

%% set the range of the travel time desired to plot
[fig1] = diffraction_image(dataseq1, T, endnoise, mySRPair,seqInitiate,neborplot,[endnoise*dt np*dt]*10^6);

disp('Click on the figure to get the travel time range')
[~,trange] = ginput(2);
close(fig1)
tinmin = floor(min(trange)*1e-6/dt+0.5);
tinmax = floor(max(trange)*1e-6/dt+0.5);
disp('Travel time range set')

%% Determine whether to plot the arrival time if it is saved previously

arrival_input = zeros(1,2);

if ~isempty(varargin) && ~isempty(varargin{1})
    arrival_input = varargin{1};
end

%% Set the reference sequence or the reference averaged signal
% by default, we take the first sequence as the reference sequence
basewiggle = squeeze(dataseq1(1,tinmin:tinmax,mySRPair.SRmap(1,1),mySRPair.SRmap(1,2)));
ref_info = []; % intialize ref_info, this is needed for the diffraction_images function
if narg>=2
    if ~isempty(varargin{2})
        ref_info = varargin{2};
        nref = length(ref_info);
        if nref == 1
            basewiggle = squeeze(dataseq1(ref_info,tinmin:tinmax,mySRPair.SRmap(1,1),mySRPair.SRmap(1,2)));
        else 
            basewiggle = mean(squeeze(dataseq1(ref_info,tinmin:tinmax,mySRPair.SRmap(1,1),mySRPair.SRmap(1,2))));
        end
    end
end

%% 
arrivalin = zeros(2,2);

%% Plot and pick
switch plotoption
    case 'global' % plot all the loaded sequences
        allwiggle = squeeze(dataseq1(1:end,tinmin:tinmax,mySRPair.SRmap(1,1),mySRPair.SRmap(1,2)));
        
        if neborplot == 1
            dataSubswiggle = [zeros(1,size(allwiggle,2)); allwiggle(2:end,1:end)-allwiggle(1:end-1,1:end)];
        else
            [dataSubswiggle] = base_signal_substraction(allwiggle,basewiggle);
        end

        % using a high-pass filter here to remove the elasto-dynamic noise
        dataSubswiggle=highpass(dataSubswiggle,0.6*1e6,Fs);
        
        figure
        imagesc((1:size(allwiggle,1))',1e6*T(tinmin:tinmax)',(dataSubswiggle(1:end,1:end))',[min(min(dataSubswiggle)),max(max(dataSubswiggle))]);
        wiggle(1e6*T(tinmin:tinmax)',(1:size(allwiggle,1))',(dataSubswiggle(1:end,1:end))','+');
        colormap('gray')
        axis ij
        hold on
        axis tight
        hold on
        plot(arrival_input(1:end,1),arrival_input(1:end,2),'k')
        hold on
        line([seqInitiate seqInitiate], [1e6*T(endnoise) 1e6*T(end)],'Color', 'red','LineStyle', '--')
        set(gcf,'Position',[100 100 1000 1000])
        ylim([min(trange) max(trange)])
        xlim([1 size(allwiggle,1)]) % this is the local sequence number
        hold on
        
        set(gcf,'NumberTitle','off');
        set(gcf,'Name','Pick arrival time for first reflection sequence');
        disp('Pick the first and the second possible sequences and arrival time');

        [in_i,arrivalin_i] = ginput(2);
        arrivalin(1:2,1) = floor(in_i+0.5);% sequence number
        arrivalin(1:2,2) = arrivalin_i;% arrival time
        
        arrival_info = arrivalin;
        plot(arrivalin(1:2,1),arrivalin(1:2,2),'ro')
        disp('End of the arrival-picking');
    case 'semi-global' % plot all the sequences waiting to be picked
        % Set the picked sequences range
        [fig2] = diffraction_image(dataseq1, T, endnoise, mySRPair,seqInitiate,neborplot,[min(trange) max(trange)],ref_info);
        disp('Click on the figure to get the range of the plotted sequence')
        [seqrange,~] = ginput(2);
        disp('Set already the range of sequence to pick')
        xin = floor(seqrange+0.5);
        xinmin = min(xin);
        xinmax = max(xin);

        close(fig2);
        
        allwiggle = squeeze(dataseq1(xinmin:xinmax,tinmin:tinmax,mySRPair.SRmap(1,1),mySRPair.SRmap(1,2)));
        if neborplot == 1
            if xinmin==1
            dataSubswiggle = [zeros(1,size(allwiggle,2)); allwiggle(2:end,1:end)-allwiggle(1:end-1,1:end)];
            else
            firstseq = squeeze(dataseq1(xinmin-1,tinmin:tinmax,mySRPair.SRmap(1,1),mySRPair.SRmap(1,2)))-allwiggle(1,1:end);
            dataSubswiggle = [firstseq; allwiggle(2:end,1:end)-allwiggle(1:end-1,1:end)];
            end
        else
            [dataSubswiggle] = base_signal_substraction(allwiggle,basewiggle);
        end
        
        % using a high-pass filter here to remove the elasto-dynamic noise
        dataSubswiggle=highpass(dataSubswiggle,0.6*1e6,Fs);
        
        figure
        imagesc((xinmin:xinmax)',1e6*T(tinmin:tinmax)',(dataSubswiggle(1:end,1:end))',[min(min(dataSubswiggle)),max(max(dataSubswiggle))]);
        wiggle(1e6*T(tinmin:tinmax)',(xinmin:xinmax)',(dataSubswiggle(1:end,1:end))','+');
        colormap('gray')
        axis ij
        hold on
        axis tight
        hold on
        plot(arrival_input(1:end,1),arrival_input(1:end,2),'k')
        set(gcf,'Position',[70 100 1000 1000])
        hold on
        plot(arrival_input(1:end,1),arrival_input(1:end,2),'k')
        hold on
        line([seqInitiate seqInitiate], [1e6*T(endnoise) 1e6*T(end)],'Color', 'red','LineStyle', '--')
        xlim([xinmin-0.5 xinmax+0.5])
        ylim([min(trange) max(trange)])
        hold on
        
        set(gcf,'NumberTitle','off');
        set(gcf,'Name','Pick arrival time for first reflection sequence');

        [in_i,arrivalin_i]=ginput(2);
        arrivalin(1:2,1) = floor(in_i+0.5);% sequence number
        arrivalin(1:2,2) = arrivalin_i;% arrival time
        
        arrival_info = arrivalin;
        plot(arrivalin(1:2,1),arrivalin(1:2,2),'ro')
        disp('End of the arrival-picking');
    otherwise
        disp('No-picking and only plotting');
        close all% close the current figures
        
        allwiggle = squeeze(dataseq1(1:end,tinmin:tinmax,mySRPair.SRmap(1,1),mySRPair.SRmap(1,2)));
        if neborplot == 1
            dataSubswiggle = [zeros(1,size(allwiggle,2)); allwiggle(2:end,1:end)-allwiggle(1:end-1,1:end)];           
        else
            [dataSubswiggle] = base_signal_substraction(allwiggle,basewiggle);
        end
        
        % using a high-pass filter here to remove the elasto-dynamic noise
        dataSubswiggle=highpass(dataSubswiggle,0.6*1e6,Fs);
        
        fig=figure
        imagesc((1:size(allwiggle,1))',1e6*T(tinmin:tinmax)',(dataSubswiggle(1:end,1:end))',[min(min(dataSubswiggle)),max(max(dataSubswiggle))]);
        if narg ==3 && ~isempty(varargin{3})
            wiggle(1e6*T(tinmin:tinmax)',(1:size(allwiggle,1))',(dataSubswiggle(1:end,1:end))','+');
        end
        colormap('gray')
        axis ij
        hold on
        axis tight
        hold on
        plot(arrival_input(1:end,1),arrival_input(1:end,2),'ko')
        hold on
        line([seqInitiate seqInitiate], [1e6*T(endnoise) 1e6*T(end)],'Color', 'red','LineStyle', '--')
        set(gcf,'Position',[100 100 1000 1000])
        ylim([min(trange) max(trange)])
        xlim([1 size(allwiggle,1)]) % this is the local sequence
        arrival_info=fig;% return the figure object
end



end
