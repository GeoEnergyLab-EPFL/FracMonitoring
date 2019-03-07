%function [Vp,Vs]=velocity_measurements()
% TODO -> write code - First AS a script and then define the proper API
%
% options- isotropy, TI, etc.

% steps
% - load data
% - pick
% -  
% output - flag for heterogeneity..... 

 
%% Function to measure sound velocity by reading two signals source-receiver at the begining

function [soundVelocity] = velocity_measurements(filepath)

% just call fuction like velocity_measurements('../data/19-01-29/133723.bin')

%% TODO:
% rewrite in function
% calc for all pairs
% average for each direction
% make soundVelocity - global variable, just function call
% a lot..

%% Length between source and receiver
   
length = 250; % normally from SampleBlock -- TBD

    
%% load data in dataSet, 16 receivers

%dataSet = cell(0,0) % need to be updated for cell
for n = 1:16 % for now only 1 direction - top/bottom -- TBD
    dataSet{n} = load_data(filepath,1,1,1); %data path
end


%take time step (dt) for data
[~, ~, jsonhdr] = load_header('../data/19-01-29/133723.json'); % read from header
Fs = jsonhdr.ActiveAcousticInfos.SamplingFrequency_MHz_*1E6; %freq of measurements 
dt = 1/Fs;  % time step

%create a timelist for easy plotting
dataTime(:,1) = 0;
for i = 1:7999
    dataTime(:,i+1) = dataTime(:,i)+dt;
end

%% reference signal read (without a sample), p-wave

dataRef = load_data('../data/19-01-11/163814.bin',1,1,1)

% figure
% plot(dataTime(:,:),abs(dataRef))


%% Check plot for dataSet
% figure
% plot(dataTime(:,:),dataSet{1})


%% Finf t_propagation with xcorr dataSet and dataRef

[acor, lag] = xcorr(dataSet{1}, dataRef);

% Figure check
% figure % plot xcorr to lag
% plot (lag*dt,abs(acor))

%t_propagation
[~,i] = max(abs(acor)); % find max of correlation, i = position
lagDiff = lag(i); [~,x]=size(lag); x=(x+1)/2; 
t_propagation = lag(i)*dt; % lagTime = timepropagation

%% Velocity mesasurment

length = 250

soundVelocity = length/(1000*t_propagation) % m/s

end
    
    
    
    
    
    
% % % % % % % % % % % % 

 % DO NOT LOOK BELOW %

% % % % % % % % % % % % 

%% 1 solution = lag = time arrival (signal, ref, [t_min, t_max]) - take global max in the time sequence



%% 2 solution = take first max>(Global max\3) of source and receiver

% 
% % \\\ example data - lvm_import for now!
% a = [1,-0.99];
% b = [1,-1];
% 
% figure
% % plot(ans.Segment1.data(:,1),ans.Segment1.data(:,2)) %source signal
% %hold on
% plot(ans.Segment1.data(:,1)*10^6,abs(ans.Segment1.data(:,3)*10^3)) %receiver signal
% %plot(dataTest.Segment1.data(:,1),dataTest.Segment1.data(:,3)) %receiver signal #1
% xlabel('Time, \mus')
% ylabel('Amplitude, mV')
% ylim([-1 1]*5) %axis tight
% set(0,'DefaultAxesColor',[1 1 1]) 
% set(0,'DefaultLineLineWidth',3) %2
% set(0,'DefaultAxesFontSize',24) %18
% set(0,'DefaultTextFontSize',22) %14
% set(0,'DefaultAxesFontWeight','bold')
% set(0,'DefaultTextFontWeight','bold')
% set(0,'DefaultLineColor','k')
% legend('Sample 1, direction 1')
% 
% 
% % \\\ find global abs max for source and receiver
% 
% %maxSource = max(abs(ans.Segment1.data(:,2)))
% %maxReceiver = max(abs(ans.Segment1.data(:,3)))
% %maxSTime = find (abs(ans.Segment1.data(:,2)) == maxSource,1,'first')
% %maxRTime = ans.Segment1.data(X,:)
% %TimeDifference = abs(abs(ans.Segment1.data(maxSTime,1)) - abs(ans.Segment1.data(maxRTime,1)))
% 
% maxGlobalRow = ans.Segment1.data(abs(ans.Segment1.data(:,3)) == max(abs(ans.Segment1.data(:,3))),:)
% maxGlobalReceiver = maxGlobalRow(1,3)
% maxGlobalRTime = maxGlobalRow(1,1)
% 
% % \\\ find first max > (Global max\3)
% 
% maxRow = ans.Segment1.data(find(abs(ans.Segment1.data(:,3)) > max(abs(ans.Segment1.data(:,3)))/3),:) % CHECK! if >maxGlobalReceiver - does not work!
% maxReceiver = maxRow(1,3)
% maxRTime = maxRow(1,1)
% 
% % looks good, repeat with source the same! 
% 
%  
%     
% %% 3 solution = pick on plot
% 
%     jj = 1; % source-receiver pair
%     figure
%     disp(['plotting source-receiver #' num2str(jj) ' amplitude over time'])
%     plot(T*1E6,dataInit2(:,:,jj,jj),T*1E6,dataInit3(:,:,jj,jj))
%     xlabel('Time (\mus)')
%     ylabel('Amplitude (a.u.)')
%     title(['source-receiver #' num2str(jj)])
%     xlim([1,2])
% 
%     [x,y] = ginput(2)
%     t_propagation = abs(x(2) - x(1))
%     soundVelocity = length/(1000*t_propagation) % m/s
% 
% end
