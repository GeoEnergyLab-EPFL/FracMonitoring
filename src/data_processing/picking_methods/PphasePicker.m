function [loc, snr_db] = PphasePicker(x, dt, type, pflag, Tn, xi, nbins, o)
%   AN AUTOMATIC P-PHASE ARRIVAL TIME PICKER
%
%   Computes P-phase arrival time in windowed digital single-component
%   acceleration or broadband velocity record without requiring threshold
%   settings. Returns P-phase arrival time in second, and signal-to-noise
%   ratio in decibel. Input waveform must be an evenly spaced vector.
%   PLEASE READ IMPORTANT NOTES BELOW
%
%   Syntax:
%      [loc, snr_db] = PphasePicker(x, dt, type, pflag, Tn, xi, nbins, o)
%
%   Input (required):
%            x = raw broadband velocity or acceleration data in
%                single-column format
%           dt = sampling interval in second (e.g., 0.005)
%         type = 'sm' or 'SM' for acceleration waveform (default bandwidth
%                = 0.1-20 Hz)
%                'wm' or 'WM' for velocity waveform (default bandwidth 7-90
%                Hz)
%                'na' or 'NA' no bandpass filtering
%
%        pflag = 'Y' or 'y' for plotting waveform with P-wave arrival time
%                marked
%                'N' or 'n' for no plot
%
%   Input (optional):
%           Tn = undamped natural period in second (default is 0.01 for
%           records sampled with 100 samples-per-second or larger; for
%           records with lower than 100 samples-per-second default
%           is 0.1 s)
%           xi = damping ratio (default is 0.6)
%        nbins = histogram bin size (default is 2/dt for
%                strong-motion acceleration and broadband velocity
%                waveforms; regional or teleseismic records may need
%                different values of bin size for better picking results)
%            o = 'to_peak' to take segment of waveform from beginning to
%                absolute peak value (recommended for fast processing)
%                'full' to take full waveform
%
%   Output:
%          loc = P-phase arrival time in second
%       snr_db = signal-to-noise ratio in decibel
%
%   Update(s):
%          2016/09/09 statelevel command is hard coded
%          2016/09/09 no bandpass filter option is added
%          2017/01/05 input normalized to prevent numerical instability
%                     as a result of very low amplitude input
%          2017/02/17 included signal-to-noise ratio computation
%
%   Example:
%          Let x be a single component strong-motion acceleration waveform
%          with 200 sample per second (dt = 0.005). The input for
%          P-phase picking will be
%
%          [loc, snr_db] = PphasePicker(x, 0.005, 'SM', 'Y');
%
%          or with optional parameters of plotting waveform (pflag =
%          'Y', Tn = 0.01, damping ratio = 0.7, histogram bin size =
%          400 and by taking segment of waveform from beginning to its
%          peak is
%
%          [loc, snr_db] = PphasePicker(x, 0.005, 'SM', 'Y', 0.01, 0.7, 400, 'to_peak')
%
%   IMPORTANT NOTE- 1: User may need to change default corner frequencies of
%   Butterworth bandpass filter as appropriate to noise level of the input
%   signal.
%
%   IMPORTANT NOTE- 2: If sampling rate of input signal is lower than 100
%   samples-per-second, use Tn = 0.1 s instead of 0.01 s to avoid numerical
%   errors. If Tn is not specified, default values will be used based on
%   the sampling rate.
%
%   Comment blocks and equation references in this function correspond to
%   the following publication:
%
%   Kalkan, E. (2016). "An automatic P-phase arrival time picker", Bull. of
%   Seismol. Soc. of Am., 106, No. 3, doi: 10.1785/0120150111
%
%   THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED
%   WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
%   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
%   NO EVENT SHALL THE COPYRIGHT OWNER BE LIABLE FOR ANY DIRECT, INDIRECT,
%   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
%   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
%   OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
%   TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
%   USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
%   DAMAGE.
%
%   Written by Dr. Erol Kalkan, P.E. (ekalkan@usgs.gov)
%   $Revision: 16.0 $  $Date: 2017/02/17 18:14:00 $
if nargin > 2
    validateattributes(x,{'double'},{'real','finite','vector'}, ...
        'PphasePicker','X');
    validateattributes(dt,{'double'},{'real','finite','scalar'}, ...
        'PphasePicker','DT');
    validatestring(type,{'sm','wm','na','SM','WM','NA'},'PphasePicker','INPUT');
else
    error('Not enough inputs.  See help documentation.');
end

if nargin > 3
    validatestring(pflag,{'Y','y','N','n'},'PphasePicker','pflag');
else
    pflag = 'n';
end

if nargin > 4
    validateattributes(Tn,{'double'},{'real','finite','scalar'}, ...
        'PphasePicker','TN');
    validateattributes(xi,{'double'},{'real','finite','scalar','<',1}, ...
        'PphasePicker','XI');
    validateattributes(nbins,{'double'},{'real','finite','scalar','>',1}, ...
        'PphasePicker','NBINS');
else
    if dt <= 0.01;
        Tn = 0.01;
        nbins = round(2/dt);
    else
        Tn = 0.1;
        nbins = 200;
    end
    xi = 0.6;
end

if nargin > 7
    validatestring(o,{'to_peak','full'},'PphasePicker','O');
else
    o = 'to_peak';
end

% User may modify the filter corner frequencies (flp and fhp) for more
% accurate picking

switch type
    % Weak-motion low- and high-pass corner frequencies in Hz
    case {'wm','WM'}
        filtflag = 1;
        %flp = 7; fhp = 90;
        flp = 0.1; fhp = 10;
        %flp = 1; fhp = 2;
        % Strong-motion low- and high-pass corner frequencies in Hz
    case {'sm','SM'}
        filtflag = 1;
        flp = 0.1; fhp = 20;
        % No bandpass filter will be applied
    case {'na','NA'}
        filtflag = 0;
        x_d = detrend(x); % detrend waveform
end
x_org = x;
% Normalize input to prevent numerical instability from very low amplitudes
x = x/max(abs(x));

% Bandpass filter and detrend waveform
if filtflag ~= 0;
    x_f = bandpass(x,flp,fhp,dt,4);
    x_d = detrend(x_f);
end

switch o
    case {'to_peak'}
        ind_peak = find(abs(x_d) == max(abs(x_d)));
        xnew = x_d(1:ind_peak);
    otherwise
        xnew = x_d;
end

% Construct a fixed-base viscously damped SDF oscillator
omegan = 2*pi/Tn;           % natural frequency in radian/second
C = 2*xi*omegan;            % viscous damping term
K = omegan^2;               % stiffness term
y(:,1) = [0;0];             % response vector

% Solve second-order ordinary differential equation of motion
A = [0 1; -K -C]; Ae = expm(A*dt); AeB = A\(Ae-eye(2))*[0;1];
for k = 2:length(xnew); y(:,k) = Ae*y(:,k-1) + AeB*xnew(k); end

veloc = (y(2,:))';          % relative velocity of mass
Edi = 2*xi*omegan*veloc.^2; % integrand of viscous damping energy

% Apply histogram method
R = statelevel(Edi,nbins);
locs = find(Edi > R(1));
indx = find(xnew(1:locs(1)-1).*xnew(2:locs(1)) < 0); % get zero crossings
TF = isempty(indx);

% Update first onset
if TF == 0;
    loc = indx(end)*dt;
else
    R = statelevel(Edi,ceil(nbins/2)); % try nbins/2
    locs = find(Edi > R(1));
    indx = find(xnew(1:locs(1)-1).*xnew(2:locs(1)) < 0); % get zero crossings
    TF = isempty(indx);
    if TF == 0; loc = indx(end)*dt; else loc = -1; end
end

% Compute SNR
if ~loc == -1;
    snr_db = -1
else
    snr_db = SNR(x,x(1:loc/dt));
end

fprintf('P-phase arrival time in second = %5.2f\n',loc);
fprintf('SNR in decibel = %5.2f\n',snr_db);

if pflag == 'Y' || pflag == 'y';
    fsz = 16; lw = 1;   % font size and line width
    figure
    set(gcf,'position',[300 283 1000 250]);
    set(gca,'TickLength',[.0025 .0025]);
    set(gca,'fontname','times','fontsize',fsz);
    T = [0:dt:(numel(x)-1)*dt];
    plot(T,x_org,'k','LineWidth',lw); grid on;
    if loc ~= -1;
        line([loc loc],[ylim],'Color','r','LineWidth',lw);
    end
    xlabel('Time, s','FontSize',[fsz+2],'fontname','times','Color','k');
    ylabel('Amplitude','FontSize',[fsz+2],'fontname','times','Color','k');
end
return

function [levels, histogram, bins] = statelevel(y,n)
ymax = max(y);
ymin = min(y)-eps;

% Compute Histogram
idx = ceil(n * (y-ymin)/(ymax-ymin));
idx = idx(idx>=1 & idx<=n);
histogram = zeros(n, 1);
for i=1:numel(idx)
    histogram(idx(i)) = histogram(idx(i)) + 1;
end

% Compute Center of Each Bin
ymin = min(y);
Ry = ymax-ymin;
dy = Ry/n;
bins = ymin + ((1:n)-0.5)*dy;

% Compute State Levels
iLowerRegion = find(histogram > 0, 1, 'first');
iUpperRegion = find(histogram > 0, 1, 'last');

iLow  = iLowerRegion(1);
iHigh = iUpperRegion(1);

% Define the lower and upper histogram regions halfway
% between the lowest and highest nonzero bins.
lLow  = iLow;
lHigh = iLow + floor((iHigh - iLow)/2);
uLow  = iLow + floor((iHigh - iLow)/2);
uHigh = iHigh;

% Upper and lower histograms
lHist = histogram(lLow:lHigh, 1);
uHist = histogram(uLow:uHigh, 1);

levels = zeros(1,2);
[~, iMax] = max(lHist(2:end));
[~, iMin] = max(uHist);
levels(1) = ymin + dy * (lLow + iMax(1) - 1.5);
levels(2) = ymin + dy * (uLow + iMin(1) - 1.5);

% Lowest histogram bin numbers for upper and lower histograms
lHist_final = (lLow + iMax(1) - 1);
uHist_final = (uLow + iMin(1) - 1);
return

function [x_f] = bandpass(c, flp, fhi, dt, n)
%  Butterworth Acausal Bandpass Filter
%
%   Syntax:
%          x_f = bandpass(c, flp, fhi, dt, n)
%
%   Input:
%            x = input time series
%          flp = low-pass corner frequency in Hz
%          fhi = high-pass corner frequency in Hz
%           dt = sampling interval in second
%            n = order
%
%   Output:
%        [x_f] = bandpass filtered signal
fnq = 1/(2*dt);              % Nyquist frequency
Wn = [flp/fnq fhi/fnq];      % Butterworth bandpass non-dimensional frequency
[b,a] = butter(n,Wn);
x_f = filtfilt(b,a,c);
return

function [r] = SNR(signal,noise)
%  Compute signal-to-noise ratio 
aps = mean(signal.^2); % average power of signal
apn = mean(noise.^2);  % average power of noise
r = 10*log10(aps/apn); % signal-to-noise ratio in decibel
return