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

tlow = 1;
tup = size(dataseq1,2)-endnoise+1;
ref = 1; 

if ~isempty(varargin) && ~isempty(varargin{1})
    ref=varargin{1};
end


if nargin>1
    if ~isempty(varargin) && ~isempty(varargin{2}) && ~isempty(varargin{3}) 
        Vp = varargin{2}(1);
        Vs = varargin{2}(2);
        
        Fs = varargin{3};
        tmid = floor((SRPairs.distances/Vp)*Fs);
        for i=1:length(SRPairs)
            if (SRPairs.wave_type(1)=='S') % when the source is a shear source
                tmid(i) = floor((SRPairs.distances/Vs)*Fs);
            end   
        end
        if ~isempty(varargin{4})
            windowsize = varargin{4};
        else
            windowsize = 100; % we set 100 as the default half window size
        end
        
        tlow = tmid-windowsize;
        tup = tmid+windowsize;   
    end
end 

if tlow<1
    tlow = 1;
end 
if tup> size(dataseq1,2)-endnoise+1
    tup=size(dataseq1,2)-endnoise+1;
end
% compute L2-norm strength
energy_all = zeros(size(dataseq1,1),size(SRPairs.SRmap,1));
for j = 1:size(SRPairs.SRmap,1)
    energy_all(:,j) = squeeze(vecnorm(dataseq1(:,tlow:tup,SRPairs.SRmap(j,1),SRPairs.SRmap(j,2)),2,2));
end
%disp(size(energy_all));
energy_all = energy_all(1:end,1:end)./energy_all(ref,1:end);

end