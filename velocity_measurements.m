
%% Function to measure sound velocity by reading two signals source-receiver at the begining


function [soundVelocity, soundVelocityPart, soundVelocityGI] = ...
    velocity_measurements(dataRefPath,dataSetPath,sequences,sources,receivers, Tmin, Tmax)

% NOW working for one-one S-R pair
% Ref signal:Swave = 162209, Pwave = 163814

% HOWTO call:
% EXAMPLE: [soundVelocity, soundVelocityPart, soundVelocityGI] = ... 
%     velocity_measurements('../data/19-01-11/163814', '../data/19-01-29/133723', 1, 1, 1, 42, 55)
% 
% (dataRefPath,dataSetPath,sequences,sources,receivers, Tmin, Tmax)
% dataRefPath -- path of the dataset, without an extantion
% dataSetPath -- path of the dataset, without an extantion
% sequences -- # of sequences in the dataSet
% sources -- # of source in the dataSet
% receivers -- # of receiver in the dataSet
% Tmin -- start of a timeframe for a direct signal, \mus
% Tmax -- end of a timeframe for a direct signal, \mus

%% TODO:
% calc for all pairs P & S waves [Vp,Vs]=velocity_measurements()
% average for each direction
% hf-experiments and Linux root

%% Length between source and receiver
   
length = 250; % normally from SampleBlock -- TBD

    
%% load data in dataSet
dataSetFilePath = strcat(dataSetPath,'.bin')
%dataSet = cell(0,0) % need to be updated for cell
%for n = 1:16 % for now only 1 pair -- all others per each side TBD
    dataSet{1} = load_data(dataSetFilePath,sequences,sources,receivers); %data path !!! N !!!
%end


%take time step (dt) for data -- normally header -- TBD
dataSetJsonPath = strcat(dataSetPath,'.json')
[~, ~, jsonhdr] = load_header(dataSetJsonPath); % read from header
Fs = jsonhdr.ActiveAcousticInfos.SamplingFrequency_MHz_*1E6; %freq of measurements 
dt = 1/Fs;  % time step

%create a timelist for easy plotting
dataTime(:,1) = 0;
for i = 1:7999
    dataTime(:,i+1) = dataTime(:,i)+dt;
end

%% reference signal read (without a sample)
dataRefFilePath = strcat(dataRefPath,'.bin')
dataRef = load_data(dataRefFilePath,1,1,1);

% figure % plot dataRef
% disp('reference signal read, source')
% plot(dataTime*1E6,dataRef,dataTime*1E6,abs(hilbert(dataRef)))


% figure % plot dataSet
% disp('plot for dataSet, receiver')
% plot(dataTime*1E6,dataSet{1})


%% Find t_propagation and Velocity with xcorr dataSet and dataRef

[acor, lag] = xcorr(dataSet{1}, dataRef);

% % Figure check
% figure % plot xcorr to lag
% disp('Xcor plot')
% plot (lag*dt,acor)

%t_propagation
[~,i] = max(acor); % find max of correlation, i = position
%lagDiff = lag(i); [~,x]=size(lag); x=(x+1)/2; 
t_propagation = lag(i)*dt; % lagTime = timepropagation

% Velocity
soundVelocity = length/(1000*t_propagation) % m/s

%% PART: Find t_propagation and Velocity with xcorr dataSetPart and dataRefPart

% set Tmin and Tmax
Tmin = Tmin*1E-6/dt; % SET Tmin for dataSet <--
Tmax = Tmax*1E-6/dt; % SET Tmax for dataSet <--

% take part of dataSet and dataRef
dataRefPart = dataRef(:,1:1000-1);
dataSetPart = dataSet{1,1}; 
dataSetPart = dataSetPart(:,round(Tmin):round(Tmax));


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
[acorPart, lagPart] = xcorr(dataSetPart, dataRefPart);

% Figure check
% figure % plot xcorr to lag
% disp('Xcor plot, PART')
% plot (lagPart*dt,acorPart)

%t_propagation
[~,j] = max(acorPart); % find max of correlation, i = position
%lagDiff = lag(i); [~,x]=size(lag); x=(x+1)/2; 
t_propagationPart = (Tmin+lagPart(j))*dt; % lagTime = timepropagation

% VelocityPart
soundVelocityPart = length/(1000*(t_propagationPart)) % m/s


%% Check with GI input

figure(gcf)
disp('GI input')
hold on
plot(dataTime*1E6,dataSet{1}*10^4) %dataSet
plot(dataTime*1E6,dataRef*10^3) %dataRef
plot(dataTime(:,round(Tmin):round(Tmax))*1E6,dataSetPart*10^4,'.g') %dataPart
% xlim([0,100]) %100/160 by X <---
xlabel('Time, \mus')
ylabel('Amplitude, mV')
legend('dataSet x10 Scale!','dataRef', 'dataSetPart for xcorr')
title('dataSet with selected Part + dataRef')


[x,y] = ginput(2)
t_propagationGI = abs(x(2) - x(1))
soundVelocityGI = length*1000/t_propagationGI % m / timeScale

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
