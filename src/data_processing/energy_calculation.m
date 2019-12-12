function energy_all  = energy_calculation(dataseq1,endnoise,SRPairs,varargin)
% Dong Liu -- 17/09/2019
% function computing energy change compared with the reference signal
% option for SR Pair: T-B all facing pairs, or one source on the top, all the
% receivers on the bottom or inverse
% option1 for reference sequence: the first local sequence by defaut
% option2 for part of the signal, Vp,Vs velocity, TODO: should add also the Vs
% velocity [Vp Vs]
% option3: the sampling frequency
% option4: the window size you take
% when calling this function, the four arguments are necessary, though
% option1 and option 4 can be [].
nargin = length(varargin);
tlow = ones(SRPairs.n_pairs,1);
tup = size(dataseq1,2).*ones(SRPairs.n_pairs,1);
ref = 1; 

if (~isempty(varargin)) && (~isempty(varargin{1}))
    ref=varargin{1};
end

if nargin>1
    if ~isempty(varargin) && ~isempty(varargin{2}) && ~isempty(varargin{3}) 

        Vp = varargin{2}(1);
        Vs = varargin{2}(2);
        Fs = varargin{3};
        tmid = floor((SRPairs.distances/Vp)*Fs);
        disp(SRPairs.wave_type(:,1));

        for i=1:SRPairs.n_pairs
            if (SRPairs.wave_type(i,1)=='S') % when the source is a shear transducer, we are more interested by its shear component
                tmid(i) = floor((SRPairs.distances(i)/Vs)*Fs);
                disp(tmid(i));
            end   
        end
        if ~isempty(varargin{4})
            windowsize = varargin{4};
        else
            windowsize = 350; % we set 350 as the default half window size
        end
        
        tlow = tmid-windowsize;
        tup = tmid+windowsize;   
    end
end 

%disp(tmid);
% build the dataseq2 for the pairs you are interested in.
% build pair by pair since the arrival time for different transducers is
% different from each other

nq = size(dataseq1,1); % number of sequences
np = size(dataseq1,2); % number of acquisition points in one sequence
dataseq2=zeros(nq,np,SRPairs.n_pairs);

% set hanning window parameter
npHanning = 50;
winHanning = hanning(npHanning*2)';
% set the filtering parameter
Fn = 0.5*Fs;
flowpass = 2E6; % cut at 2 MHz
[b, a] = butter(3,flowpass/Fn);

for i_d=1:SRPairs.n_pairs
    localseq = dataseq1(:,:,SRPairs.SRmap(i_d,1),SRPairs.SRmap(i_d,2));
    Twindow = [zeros(1,tlow(i_d)-npHanning) winHanning(1:npHanning) ones(1,tup(i_d)-tlow(i_d))...
    winHanning(npHanning+1:end) zeros(1,np-(tup(i_d)+npHanning))]';
    datasetWin = localseq.*Twindow';
    dataLowFilt = filtfilt(b,a,datasetWin');
    dataseq2(:,:,i_d)= dataLowFilt';
end


% Twindow = [zeros(1,tlow-npHanning) winHanning(1:npHanning) ones(1,tup-tlow)...
%     winHanning(npHanning+1:end) zeros(1,np-(tup+npHanning))]';
% 
% % window signal with it
% Twindow3D = permute(repmat(Twindow,1,nq,npairs),[2,1,3]);% change the dimension
% disp(num2str(size(Twindow3D)));
% disp(num2str(size(dataseq1)));
% datasetWin = dataseq1.*Twindow3D;
% 
% 
% datasetWin = permute(datasetWin,[2 1 3 4]);% we need put the time-axis on the first dimension to make filtfilt work.
% dataLowFilt = filtfilt(b,a,datasetWin);
% dataseq2 = permute(dataLowFilt,[2 1 3 4]);

% compute L2-norm strength
energy_all = zeros(size(dataseq2,1),size(SRPairs.SRmap,1));
for j = 1:size(SRPairs.SRmap,1)
    energy_all(:,j) = squeeze(vecnorm(dataseq2(:,tlow:tup,j),2,2));
end
%disp(size(energy_all));
energy_all = energy_all(1:end,1:end)./energy_all(ref,1:end);

end