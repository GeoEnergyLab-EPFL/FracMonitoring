function [arrival_info]=diffractions_picking(activeInfo,dataseq1, T, endnoise, mySRPair, seqInitiate, plotoption, varargin)
% Dong Liu -- 01/10/2019
% the detailed plot for picking the arrival time for a certain pair
% in order to have a closer look
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
% plotoption: one only has three choices for plotoptions:
% global,semi-global and local
% optional argument1: give the (nx2) matrix (previously saved data for this pair), in
% case of the first-time picking, this matrix value is all 0--> zeros(1,2)
% optional argument2: give the reference sequence number 
% or the sequence ranges where the the averaged reference
% sequence is calculated from
% optional argument3: indicate whether to plot the wiggle plot or not for
% the plotting part
%% get the active acoustic info
np = activeInfo.NumberOfPoints;
Fs = activeInfo.SamplingFrequency_MHz_*1E6;
dt = 1/Fs;  % time step

%% set the range of the travel time desired to plot
[fig1] = diffraction_image(dataseq1, T, endnoise, mySRPair,seqInitiate,0,[endnoise*dt np*dt]*10^6);

disp('Click on the figure to get the travel time range')
[~,trange] = ginput(2);
close(fig1)
tinmin = floor(min(trange)*1e-6/dt+0.5);
tinmax = floor(max(trange)*1e-6/dt+0.5);
disp('Travel time range set')

%% Determine whether to plot the arrival time if it is saved previously

narg = length(varargin);
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

%% Set the picked sequences range
[fig2] = diffraction_image(dataseq1, T, endnoise, mySRPair,seqInitiate,0,[min(trange) max(trange)],ref_info);
disp('Click on the figure to get the range of the plotted sequence')
[seqrange,~] = ginput(2);
disp('Set already the range of sequence to pick')
xin = floor(seqrange+0.5);
xinmin = min(xin);
xinmax = max(xin);

close all;

arrivalin = zeros(xinmax-xinmin+1,2);

%% Plot and pick
switch plotoption
    case 'global' % plot all the loaded sequences
        allwiggle = squeeze(dataseq1(1:end,tinmin:tinmax,mySRPair.SRmap(1,1),mySRPair.SRmap(1,2)));
        [dataSubswiggle] = base_signal_substraction(allwiggle,basewiggle);
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
        xticks(1:5:size(allwiggle,1)) % this is true only when we read the sequences continously from 1
        hold on
        for in_i = xinmin:xinmax
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',['Pick arrival time for Seq' num2str(in_i)]);
            if in_i > xinmin
                plot(arrivalin(1:in_i-xinmin,1),arrivalin(1:in_i-xinmin,2),'ro')
            end
            [~,arrivalin_i] = ginput(1);
            arrivalin(in_i-xinmin+1,1) = in_i;% sequence number
            arrivalin(in_i-xinmin+1,2) = arrivalin_i;% arrival time
        end
        arrival_info = arrivalin;
        plot(arrivalin(1:end,1),arrivalin(1:end,2),'r')
        disp('End of the arrival-picking');
    case 'semi-global' % plot all the sequences waiting to be picked
        allwiggle = squeeze(dataseq1(xinmin:xinmax,tinmin:tinmax,mySRPair.SRmap(1,1),mySRPair.SRmap(1,2)));
        [dataSubswiggle] = base_signal_substraction(allwiggle,basewiggle);
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
        xticks(1:5:size(allwiggle,1)) % this is true only when we read the sequences continously from 1
        hold on
        for in_i = xinmin:xinmax
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',['Pick arrival time for Seq' num2str(in_i)]);
            if in_i>xinmin
                plot(arrivalin(1:in_i-xinmin,1),arrivalin(1:in_i-xinmin,2),'ro')
            end
            [~,arrivalin_i]=ginput(1);
            arrivalin(in_i-xinmin+1,1) = in_i;% sequence number
            arrivalin(in_i-xinmin+1,2) = arrivalin_i;% arrival time
        end
        arrival_info = arrivalin;
        plot(arrivalin(1:end,1),arrivalin(1:end,2),'r')
        disp('End of the arrival-picking');
    case 'local' % plot one sequence and its neighboring sequences(one set 10 here), replot the figure each time
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        nbspare=10; % sequences shown in the wiggle before and after the picking sequence, this value can be changed if needed
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        for in_i = xinmin:xinmax
            allwiggle = squeeze(dataseq1(in_i-nbspare:in_i+nbspare,tinmin:tinmax,mySRPair.SRmap(1,1),mySRPair.SRmap(1,2)));
            [dataSubswiggle] = base_signal_substraction(allwiggle,basewiggle);

            fig_ini = in_i;
            figure(fig_ini)
            hold on
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',['Pick arrival time for Seq' num2str(in_i)]);
            imagesc((in_i-nbspare:in_i+nbspare)',1e6*T(tinmin:tinmax)',(dataSubswiggle(1:end,1:end))',[min(min(dataSubswiggle)),max(max(dataSubswiggle))]);
            colormap('gray')
            axis ij
            hold on
            axis tight
            wiggle(1e6*T(tinmin:tinmax)',(in_i-nbspare:in_i+nbspare)',(dataSubswiggle(1:end,1:end))','+');
            hold on
            if in_i > xinmin
            plot(arrivalin(1:in_i-xinmin,1),arrivalin(1:in_i-xinmin,2),'ro')
            end
            hold on
            line([seqInitiate seqInitiate], [1e6*T(endnoise) 1e6*T(end)],'Color', 'red','LineStyle', '--')
            xlim([in_i-nbspare in_i+nbspare])
            ylim([min(trange) max(trange)])
            [~,arrivalin_i] = ginput(1);
            arrivalin(in_i-xinmin+1,1) = in_i;% sequence number
            arrivalin(in_i-xinmin+1,2) = arrivalin_i;% arrival time
            close(fig_ini)
        end
        arrival_info = arrivalin;
        disp('End of the arrival-picking');
        close all% close the current figures
        [~] = diffraction_image(dataseq1, T, endnoise, mySRPair,seqInitiate,0,[min(trange) max(trange)],ref_info);
        hold on 
        plot(arrivalin(1:end,1),arrivalin(1:end,2),'r')
    otherwise
        disp('No-picking and only plotting');
        close all% close the current figures
        
        allwiggle = squeeze(dataseq1(1:end,tinmin:tinmax,mySRPair.SRmap(1,1),mySRPair.SRmap(1,2)));
        [dataSubswiggle] = base_signal_substraction(allwiggle,basewiggle);
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
        plot(arrival_input(1:end,1),arrival_input(1:end,2),'k')
        hold on
        line([seqInitiate seqInitiate], [1e6*T(endnoise) 1e6*T(end)],'Color', 'red','LineStyle', '--')
        set(gcf,'Position',[100 100 1000 1000])
        ylim([min(trange) max(trange)])
        xticks(1:5:size(allwiggle,1)) % this is true only when we read the sequences continously from 1
        hold on
        plot(arrival_input(1:end,1),arrival_input(1:end,2),'k')
        arrival_info=fig;% return the figure object
end



end
