%% Function to measure sound velocity by reading two signals source-receiver at the begining


function [VmatrixOut] = ...
    velocity_advanced(VmatrixIn, TimeP, TimeS, method)

% Ref signal:Swave = 162209, Pwave = 163814

% HOWTO call: [VmatrixOut] = velocity_advanced(VmatrixIn, {30,50}, {60,80}, ["Full", "Part", "GI", "AIC"]); % Vmatrix, {Pwave Tmin, Pwave Tmax}, {S wave Tmin, Swave Tmax}, ["Full", "Part", "GI", "AIC", "STA_LTA"], dt

% VmatrixIn -- input parameters from data before
% TimeP -- range of time Tmin and Tmax for P wave \mus
% TimeS -- range of time Tmin and Tmax for S wave \mus
% method -- ["Full", "Part", "GI", "AIC"], which method do you want to start -
% Full data correlation, selected Part, Graphical Input, AIC

%Plot setting
set(0,'DefaultAxesColor',[1 1 1]) 
set(0,'DefaultLineLineWidth',3) %2
set(0,'DefaultAxesFontSize',24) %18
set(0,'DefaultTextFontSize',22) %14
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultTextFontWeight','bold')
set(0,'DefaultLineColor','k')

% ADDPATH
addpath('../src/data_processing/');
addpath('../src/data_load/');
addpath('../src/utilities/');
addpath('../src/properties/');
addpath('../src/inversion/');
addpath('../src/forward/');
addpath('../src/experiment_configuration/');
addpath('../src/data_processing/picking_methods/');


% Time matrix for graphs
dataTime(:,1) = 0;
for i = 1:7999
    dataTime(:,i+1) = dataTime(:,i)+VmatrixIn.dt;
end

        
% Plot PP and SS to figure out time frames
figure %(gcf)
disp('Take a look on PP and SS plots')
hold on
plot(dataTime*1E6,VmatrixIn.dataset(:,find(VmatrixIn.Ttype == "PP",1))*10^3) %dataSet PP
plot(dataTime*1E6,VmatrixIn.dataset(:,find(VmatrixIn.Ttype == "SS",1))*10^4) %dataSet SS
% xlim([0,100]) %100/160 by X <---
xlabel('Time, \mus')
ylabel('Amplitude, mV')
title('PP and SS plots')
legend('PP','SS x10')
set(gca,'box','on')

disp('Click on the plot to go next:');
[~]=ginput(1);
% close
       
 % figure % plot dataRef
 % disp('reference signal read, source')
 % plot(dataTime*1E6,dataRef)

% [~,nr] = size (dataSet(1,:)); %number of receivers



%% Method 1: Full xcorr
% Find t_propagation and Velocity with xcorr dataSet and dataRef

if (find(contains(method, 'Full')) > 0) % if method is selected
    for ii = 1:VmatrixIn.N_pairs % go over all pairs
        if strcmp("PP",VmatrixIn.Ttype(ii)) | strcmp("PS",VmatrixIn.Ttype(ii)) % if PP or PS transducers - check with P ref signal, else with S

            [acor, lag] = xcorr(VmatrixIn.dataset(:,ii), VmatrixIn.dataRef_P); % dataRef depends on the Transducer type -- for 'PP' or 'PS'
            
            [~,i] = max(acor); % find max of correlation, i = position
            t_propagation = lag(i)*VmatrixIn.dt; % lagTime = time propagation

            % Velocity P
            VmatrixOut.VpFull(ii,1) = VmatrixIn.Length(ii) / t_propagation; % m/s
        
        else % is SS or SP types
            [acor, lag] = xcorr(VmatrixIn.dataset(:,ii), VmatrixIn.dataRef_S); % for 'SS' and 'SP'
            
            [~,i] = max(acor); % find max of correlation, i = position
            t_propagation = lag(i)*VmatrixIn.dt; % lagTime = time propagation

            % Velocity S
            VmatrixOut.VsFull(ii,1) = VmatrixIn.Length(ii) / t_propagation % m/s
        end
        % % Figure check xcorr to lag
        % figure
        % disp('Xcor plot')
        % plot (lag*VmatrixIn.dt,acor)

    end
    
else % if Full method is not selected
        VmatrixOut.VpFull = 'Not selected';
        VmatrixOut.VsFull = 'Not selected';

end

VmatrixOut.VpFull_median = median(nonzeros(VmatrixOut.VpFull)); % median data over VpFull
VmatrixOut.VsFull_median = median(nonzeros(VmatrixOut.VsFull)); % median data over VsFull

%% Method 2: PART xcorr 
% Find t_propagation and Velocity with xcorr dataSetPart and dataRefPart

if (find(contains(method, 'Part')) > 0)

        % set Tmin and Tmax
        TimeP{1} = TimeP{1}*1E-6/VmatrixIn.dt; % SET Tmin for P dataSet <--
        TimeP{2} = TimeP{2}*1E-6/VmatrixIn.dt; % SET Tmax for P dataSet <--
        TimeS{1} = TimeS{1}*1E-6/VmatrixIn.dt; % SET Tmin for S dataSet <--
        TimeS{2} = TimeS{2}*1E-6/VmatrixIn.dt; % SET Tmax for S dataSet <--
    
        % take part of dataRef for S and P
        VmatrixIn.dataRef_SPart = VmatrixIn.dataRef_S(1:1000-1,:);
        VmatrixIn.dataRef_PPart = VmatrixIn.dataRef_P(1:1000-1,:);
        

        for ii = 1:VmatrixIn.N_pairs
    
            if strcmp("PP",VmatrixIn.Ttype(ii)) | strcmp("PS",VmatrixIn.Ttype(ii)) % if PP or PS Ttype
                [acorPart, lagPart] = xcorr(VmatrixIn.dataset(round(TimeP{1}):round(TimeP{2}),ii), VmatrixIn.dataRef_PPart); % dataRef depends on the Transducer type -- for 'PP' or 'PS'

            [~,j] = max(acorPart); % find max of correlation, i = position
            %lagDiff = lag(i); [~,x]=size(lag); x=(x+1)/2; 
            t_propagationPart = (TimeP{1}+lagPart(j))*VmatrixIn.dt; % lagTime = timepropagation

            % VelocityPart P
            VmatrixOut.VpPart(ii,1) = VmatrixIn.Length(ii) / t_propagationPart; % m/s

            else % if SS or SP Ttype
                [acorPart, lagPart] = xcorr(VmatrixIn.dataset(round(TimeS{1}):round(TimeS{2}),ii), VmatrixIn.dataRef_SPart); % for 'SS' and 'SP'

            [~,j] = max(acorPart); % find max of correlation, i = position
            %lagDiff = lag(i); [~,x]=size(lag); x=(x+1)/2; 
            t_propagationPart = (TimeS{1}+lagPart(j))*VmatrixIn.dt; % lagTime = timepropagation

            % VelocityPart S
            VmatrixOut.VsPart(ii,1) = VmatrixIn.Length(ii) / t_propagationPart; % m/s

            end 
            
        end
else
    VmatrixOut.VpPart = 'Not selected';
    VmatrixOut.VsPart = 'Not selected';

end

VmatrixOut.VpPart_median = median(nonzeros(VmatrixOut.VpPart));
VmatrixOut.VsPart_median = median(nonzeros(VmatrixOut.VsPart));


%% Method 3: GI
%Check with GI input

if (find(contains(method, 'GI')) > 0)
    
    for ii = 100:101 %VmatrixIn.N_pairs >>>NOT TO DIE PICKING ALL -- NOW JUST 2 PAIRS<<<
% % 100 - SS, 101 - PP

        
        nowPlot = figure %(gcf)
        disp('GI input')
        xlim([0,100]) %100/160 by X <---
        xlabel('Time, \mus')
        ylabel('Amplitude, mV')
        set(gca,'box','on')
        % title('Amlitude of Source and Receiver')
        hold on
        
        plot(dataTime*1E6,VmatrixIn.dataset(:,ii)*10^3) %dataSet
        
        if strcmp("PP",VmatrixIn.Ttype(ii)) | strcmp("PS",VmatrixIn.Ttype(ii)) % >>> can be improved -- first to plot and after check fot Ttype, devide into different coloumns <<<
            plot(dataTime*1E6,VmatrixIn.dataRef_P*10^3); %dataRef_P
            plot(dataTime(:,round(TimeP{1}):round(TimeP{2}))*1E6,VmatrixIn.dataset(round(TimeP{1}):round(TimeP{2}),ii)*10^3,'.g') %dataPart P
            legend('Initial signal','Reference signal', 'Selected data')

            figure(nowPlot)
            [x,y] = ginput(2)
            t_propagationGI = abs(x(2) - x(1))
            VmatrixOut.VpGI(ii,1) = VmatrixIn.Length(ii) * 10^6 / t_propagationGI % m / timeScale
            % soundVelocityGI = 0;            
        else
            plot(dataTime*1E6,VmatrixIn.dataRef_S*10^3); %dataRef_S
            plot(dataTime(:,round(TimeS{1}):round(TimeS{2}))*1E6,VmatrixIn.dataset(round(TimeS{1}):round(TimeS{2}),ii)*10^3,'.g') %dataPart S
            legend('Initial signal','Reference signal', 'Selected data')

            figure(nowPlot)
            [x,y] = ginput(2)
            t_propagationGI = abs(x(2) - x(1))
            VmatrixOut.VsGI(ii,1) = VmatrixIn.Length(ii) * 10^6 / t_propagationGI % m / timeScale
            % soundVelocityGI = 0;           
        end
        
        hold off
        close

    end
else
    VmatrixOut.VpGI = 'Not selected';
    VmatrixOut.VsGI = 'Not selected';

end

VmatrixOut.VpGI_median = median(nonzeros(VmatrixOut.VpGI));
VmatrixOut.VsGI_median = median(nonzeros(VmatrixOut.VsGI));



%% Method 4: AIC criteria 

% works only to detect P as higher amplitude? 


if (find(contains(method, 'AIC')) > 0)
    
    endnoise = 100; % burning period for EM noise from the source excitation
    % brute force automatic peakinf of first wave arrival using AIC
    % Brice April 15 2019   -> to be  improved and then move to a function

    % iterates on pair in the map
    % we use the first acquisition sequence only here - Initial velocity. 
    type=[];
    slowness=zeros(VmatrixIn.N_pairs,3);
    for ii=1:VmatrixIn.N_pairs
        signal=VmatrixIn.dataset(:,ii)

        type=[type; VmatrixIn.Ttype(ii,:)];

        % figure
        % plot(signal,'.')
        % aic ?
        tic
        aic = AIC(signal);
        toc

        % figure
        % plot(aic)
        %
        % figure
        % plot(diff(aic))


        % now let's peak the arrival using 3 different criteria
        [f k]=min(aic(endnoise:end-10));  % minimum of AIC - not reliable (with multiple arrival)    
        [f2 k2]=max(diff(aic(endnoise:end-10))); % max of AIC derivatives (>0) -> good for first arrival in most cases    
        [f3 k3 ]=max(signal(endnoise:end)); % maximum signal peak -> lower velocity 

        dtpick=[k+endnoise k2+endnoise k3+endnoise];
        slownessAIC(ii,:)=dtpick.*VmatrixIn.dt/VmatrixIn.Length(ii); % estimate slowness 1/Vel
    end
    
    VmatrixOut.VAIC=1./slownessAIC; % corresponding velocity
    % some plottings
%     figure
%     plot(VmatrixOut.soundVelocityAIC);
%    
    VmatrixOut.VAIC_median(1,1) = median(VmatrixOut.VAIC(:,1));
    VmatrixOut.VAIC_median(1,2) = median(VmatrixOut.VAIC(:,2));
    VmatrixOut.VAIC_median(1,3) = median(VmatrixOut.VAIC(:,3));

    
    
end
    

%% median of all (Full, Part, GI)

VmatrixOut.VpAverage = [VmatrixOut.VpFull_median; VmatrixOut.VpPart_median; VmatrixOut.VpGI_median] % 101 means no value
VmatrixOut.VsAverage = [VmatrixOut.VsFull_median; VmatrixOut.VsPart_median; VmatrixOut.VsGI_median] % 101 means no value



end



% DO NOT LOOK BELOW





% %% look at one particular combinaison of the source-receiver pair in the map

% ii=2;
% signal=VmatrixIn.dataset(:,ii)
% 
% aic = AIC(signal);
% 
% figure 
% title('Signal')
% plot(signal,'.')
% 
% figure
% title('AIC')
% plot(aic)



%% Method 5: STA/LTA 

% BY DONG

% function tindex=arrival_time(data1seq)

% abs_v = abs(data1seq);% it should be a 1xN array
% l_v=size(data1seq,2);% length of the datapoint
% trig_array = zeros(1,2);% initialization
% ntrig=0;
% lta_calc_flag = 0; % to test if the full calculation is needed or not?
% 
% % set all the controlling parameters
% l_sta = round(4/160*8000);     % STA window length 4 us
% l_lta = round(20/160*8000);     % LTA window length 20 us
% th_on = 1.2;        % Trigger on when sta_to_lta exceeds this theshold
% th_off = 1.;     % Trigger off when sta_to_lta drops below threshold
% min_dur = round(0/160*8000);   % Any triggers shorter than min_dur are discarded
% 
% ratio_trig = zeros(l_v-l_lta,1); % record the stalta_ratio for each i value
% 
% % i is the primary reference point (right end of STA/LTA window)
% i = l_lta+1;
% while i <= l_v % START STA_LTA MAIN LOOP
%     if (lta_calc_flag == 0)
%       lta_sum = 0;
%       sta_sum = 0;
%       for j = i-l_lta:i-1              % Loop to compute LTA & STA
%          lta_sum = lta_sum + abs_v(j); % Sum LTA window
%          if (i - j) <= l_sta           % Sum STA window (right side of LTA)
%             sta_sum = sta_sum + abs_v(j);
%          end
%       end
%       lta_calc_flag = 1; % after the first calculation, the following is just a moving window, no need to recalculate everything again
%     else
%       lta_sum = lta_sum - abs_v(i-l_lta-1) + abs_v(i-1);
%       sta_sum = sta_sum - abs_v(i-l_sta-1) + abs_v(i-1);
%     end
%     lta = lta_sum/l_lta;
%     sta = sta_sum/l_sta;
%     sta_to_lta = sta/lta;
%     
%    if (sta_to_lta > th_on)
%       j = i;   % Set secondary reference point = primary
%       g = 0;   % l_lta growth, only used if LTA growing
%       while (sta_to_lta > th_off)
%          j = j+1;
%          if j < l_v
%             sta_sum = sta_sum - abs_v(j-l_sta-1) + abs_v(j-1);
%             
% %             switch lta_mode
% %                case 'frozen'
% %                   % LTA is good just the way it is
% %                case 'continuous'
%             % Add new data point, remove oldest data point
%             % we assume here the windowsize of LTA moves together with the
%             % STA  to the right
%              lta_sum = lta_sum - abs_v(j-l_lta-1) + abs_v(j-1);
% %                case 'grow'
% %                   % Add new data point, increase
% %                   lta_sum = lta_sum + abs_v(j-1);
% %                   l_lta = l_lta + 1;
% %                   g = g+1;
% %             end
%             sta = sta_sum/l_sta;
%             lta = lta_sum/l_lta;
%             sta_to_lta = sta/lta;
%          else
%             sta_to_lta = 0; % Force trigger off (end of data)
%          end
%       end
%       duration = (j-i); % span from trigger on to trigger off
%       l_lta = l_lta-g;
%       if duration >= min_dur % If duration < min_dur then skip it
%          trig_t = i-l_sta;  % Beginning of STA window during trigger on
%          end_t  = j;        % End of STA window during trigger off
%          ntrig = ntrig + 1; % Event counter
%          trig_array(ntrig,:) = [trig_t, end_t];
%       end
%       lta_calc_flag = 0;
%    end
%    ratio_trig(i-l_lta,1)=sta_to_lta;
%    i = i + 1;
% end
% 
% if (trig_array(1,1)==0)&&(trig_array(1,2)==0)
%     disp('No events detected')
%     tindex = 0;
%     return
% end
% 
% tindex = trig_array(1,1);
% disp(tindex);
% disp(["t index is ", num2str(tindex), num2str(sta_to_lta)]);
% figure
% plot((1:size(ratio_trig,1))'/8000*160+l_lta*160/8000,ratio_trig)
% xlabel('Travel time (\mu s)')
% ylabel('STA/LTA')



%% Method 6: PAI-K picker

% BY DONG

% function [ind] = pai_kpicker(x,M,o)

% Analysis of kurtosis (pick) by windows

% the PAI-K picker can be seriously affected by a low SNR, spikes, a strong S wave and large amplitude at the tails.

%   Input:
%          x = raw data in single-column format
%          M = the window size, should be smaller than the signal length
%   Output:
%          ind = Index of pick in x; pick is the min(AIC(n)) + 1
%

% x = x - median(x); % remove median of x and window
% switch o
%     case {0,'1','to_peak'}
%         ind_peak = find(abs(x) == max(abs(x)));
%         xnew = x(1:ind_peak);
%     otherwise
%         xnew = x;
% end
% junk=zeros(size(xnew,1),1);
% for k=M:size(xnew,1)
%     junk(k,1) = kfunc(xnew,k,M);
% end
% junk(1:M-1)=min(junk);
% if junk ~= 0
%     ind = find(junk == max(junk)); % corresponds to the maximum value
% else
%     ind = 0;
% end
% 
% % show the plot of the AIC index and together with the signal
% figure
% plot(1:size(junk,1),junk)
% xamp=max(junk)-min(junk);
% hold on
% plot(1:size(xnew,1),xnew/max(xnew)*xamp+mean(junk))
% hold on
% line([ind ind], [min(junk) max(junk)],'Color', 'black','LineStyle', '--')
% 
% return
% 
% 
% 
% % function [a] = kfunc(x,k,M)
%     if ~isempty(x)
%     xnew=x(k-M+1:k);
%     mk=mean(xnew);
%     aup=0;
%     adown=0;
%     for i=k-M+1:k
%         aup=aup+(x(i)-mk)^4;
%         adown=adown+(x(i)-mk)^2;
%     end
%     a=(M-1)*aup/adown-3;
% else
%     a = 0;
% end
% return


%% Method 7: P-phase Picker

% looks like works bad with high Hz

% BY DONG

% data_PP = VmatrixIn.dataset(:,1);
% [loc, snr_db] = PphasePicker(data_PP, VmatrixIn.dt, 'na', 'Y');

% function [loc, snr_db] = PphasePicker(x, VmatrixIn.dt, type, pflag, Tn, xi, nbins, o)


%% Method 8: auto-picking

% look for kortois, freq analysis, depends on rock autopick, 


%% Method 9: ML picking



%% Another?

% ???? - Deep learnnig like in NLP - voice recognition













%% load data in dataSet
%%q dataSetFilePath = strcat(dataSetPath,'.bin')
%dataSet = cell(0,0) % need to be updated for cell
%for n = 1:16 % for now only 1 pair -- all others per each side TBD
%%q    dataSet{1} = load_data(dataSetFilePath,sequences,sources,receivers); %data path !!! N !!!
%end
%    dataSet = dataset(:,1); %data path !!! N !!!


%% reference signal read (without a sample)
%%q dataRefFilePath = strcat(dataRefPath,'.bin')
%%q dataRef = load_data(dataRefFilePath,1,1,1);

% dataRef = dataRef;


% figure % plot dataRef
% disp('reference signal read, source')
% plot(dataTime*1E6,dataRef,dataTime*1E6,abs(hilbert(dataRef)))


% figure % plot dataSet
% disp('plot for dataSet, receiver')
% plot(dataTime*1E6,dataSet)


%% choose all pairs with good signal
% 
% % Quality Control of acoustic signal: SNR matrix on all pairs.
% endnoise = 100; % burning period for EM noise from the source excitation
% % [normAll, h_SNR]  = signal_strength(dataSet),endnoise)
% 
% 
% 
% 
% % Selection of source - receiver pairs to compute velocities at t=0
% 
% % We choose all pairs from opposite platen 
% oppSR=AllPairsOppositePlattens(myTransducers,myPlattens);
% 
% fig_b=plotblockwithplattens(myBlock,myPlattens);
% %
% fig_handle=plotdirectrays(oppSR,fig_b);
% 
% 
% 
% 
% % get the S-R pairs  with sufficient snr
% 
% snrv=zeros(length(oppSR.SRmap),1); % here I do a loop because Matlab is stupid / not well designed for array access
% for i=1:length(snrv)
%     snrv(i)=normAll(oppSR.SRmap(i,1),oppSR.SRmap(i,2));
% end
% 
% figure 
% plot(snrv)
% 
% [r c ~]=find(snrv>2);  % adjust threshold of snr above which we look at the data - we are strict here
% 
% snrvf=snrv(r);
% 
% % create a new object source-receiver pairs containing the pairs from
% % opposite platen with a sufficient SNR 
% submap= oppSR.SRmap(r,:);
% subwave=oppSR.wave_type(r,:);
% 
% SRdirect=SourceReceiverPairs(myTransducers,myPlattens,submap,subwave);
% fig_b=plotblockwithplattens(myBlock,myPlattens);
% 
% fig_handle=plotdirectrays(SRdirect,fig_b);
% 
% 
% 
%%


            % PLOT dataSet + dataSetPart on 1 figure
            % figure % plot dataSet
            % plot(dataTime*1E6,dataSet{1})
            % hold on
            % plot(dataTime(:,round(Tmin):round(Tmax))*1E6,dataSetPart,'.r')
            % legend('dataSet','dataSetPart')
            % title('dataSet + dataSetPart')

            % figure % plot dataRefPart
            % disp('Xcor part, dataRefPart')
            % plot(dataTime(:,1:1000-1),dataRefPart)

            % figure % plot dataSetPart
            % plot(dataTime(:,Tmin:Tmax),dataSetPart)


            % xcorr with dataRefPart and dataSetPart

             %figure check xcorr to lag
             %figure 
             %disp('Xcor plot, PART')
             %plot (lagPart*VmatrixIn.dt,acorPart)

            %t_propagation
    %         [~,j] = max(acorPart); % find max of correlation, i = position
    %         %lagDiff = lag(i); [~,x]=size(lag); x=(x+1)/2; 
    %         t_propagationPart = (Tmin+lagPart(j))*VmatrixIn.dt; % lagTime = timepropagation
    % 
    %         % VelocityPart
    %         VmatrixOut.soundVelocityPart(ii,1) = VmatrixIn.Length(ii) / t_propagationPart % m/s
    %     end
