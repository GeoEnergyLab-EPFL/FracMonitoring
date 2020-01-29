function [figs]=diffraction_image(dataseq1, T, endnoise, mySRPair, seqInitiate, diffindicator, plotrange, varargin)
% Dong Liu -- 20/09/2019
% to plot the diffraction curve (all the loaded sequences) of a list of 
% transducer pairs in order to choose the best ones to pick the 
% diffraction arrival
%
% dataseq1: the loaded acoustic data
% mySRPair: chosen transducer pairs object
% seqInitiate: the guessed number of the sequence corresponding to the
% fracture initiation, obtained from the presssure responses
% diffindicator: bool value, if the value==1 plot the difference between the
% neighboring sequences
% plotrange: the plot range for the travel time
% optional argument 1: set the reference sequence
% optional argument 2: set whether to plot the wiggle plot
% optional argument 3: active acoustic info: [Fs F_excitation]
% optional argument 4: filtering range for the signals

%% set the reference sequence
ref_nb=1;
narg=length(varargin);

if narg>=1
    if ~isempty(varargin) && ~isempty(varargin{1})
        ref_nb=varargin{1};
    end 
end

%% determine whether to filter the signals
fil_range=[]; 
n_fil=length(fil_range);
if narg>=3 && ~isempty(varargin{3})
    activeInfo=varargin{3};
    Fs = activeInfo(1); % sampling frequency
    F_excitation= activeInfo(2); % excitation frequency
    if narg>=4 && ~isempty(varargin{4})
        fil_range=varargin{4};
        n_fil=length(fil_range);
        if n_fil>=1
            disp('A high-pass filter will be applied');
            if n_fil==2
                disp('A low-pass filter will be also applied');
            end
        end
    end
end

%% Diffraction plot for all pairs
figs=zeros(mySRPair.n_pairs,1);
for j=1:mySRPair.n_pairs
    if length(ref_nb)>1
        basesignal=mean(squeeze(dataseq1(ref_nb,endnoise:end,mySRPair.SRmap(j,1),mySRPair.SRmap(j,2))));
    else
        basesignal=squeeze(dataseq1(ref_nb,endnoise:end,mySRPair.SRmap(j,1),mySRPair.SRmap(j,2)));
    end
    
    allsignal=squeeze(dataseq1(1:end,endnoise:end,mySRPair.SRmap(j,1),mySRPair.SRmap(j,2)));
    [dataSubs] = base_signal_substraction(allsignal,basesignal);
    
    if diffindicator == 1
        newdataSubs=[zeros(1,size(dataSubs,2)); dataSubs(2:end,1:end)-dataSubs(1:end-1,1:end)];
        % newdataSubs=dataSubs-newdataSub; 
        % dataSubs=a(n)-a(0); newdataSubs=a(n)-a(n-1)
        dataSubs=newdataSubs;
    end
    
    % using a high-pass filter here to remove the elasto-dynamic noise
    if n_fil>=1
        dataSubs=highpass(dataSubs,fil_range(1),Fs);
        if n_fil==2
            dataSubs=lowpass(dataSubs,fil_range(2),Fs);
        end
    end
    
%      dataSubs=highpass(dataSubs',0.001,1);
%      dataSubs=dataSubs';
%      dataSubs=highpass(dataSubs,400*1e3,50*1e6);
%      dataSubs=lowpass(dataSubs,50000*1e3,50*1e6);

    fig_j=j;
    figure(fig_j)
    hold on;
    % !!!!! please do not change the title here, the chosen pairs will take the
    % fig number from the title
    title(['Fig' num2str(j) '  : S' num2str(mySRPair.SRmap(j,1)) '-R' ...
        num2str(mySRPair.SRmap(j,2)) ' of type '  num2str(mySRPair.wave_type(j,:))...
        ' distance='  num2str(mySRPair.distances(j))],'fontsize',14);

    imagesc(1:size(dataSubs,1),1e6*T(endnoise:end),dataSubs',[0.9*min(min(dataSubs)),0.9*max(max(dataSubs))])
    axis ij
    axis tight
    colormap('gray')
    colorbar
    hold on
    if narg>=2 && ~isempty(varargin{2})
        wiggle(1e6*T(endnoise:end)',(1:size(allsignal,1))',(dataSubs(1:end,1:end))','+');
    end
    hold on

    xlabel('Local sequence ()')
    ylabel('Travel time (\mu s)')
    hold on
    line([seqInitiate seqInitiate], [1e6*T(endnoise) 1e6*T(end)],'Color', 'red','LineStyle', '--')
    ylim(plotrange);
    
    figs(j,1)=fig_j;

end

end
