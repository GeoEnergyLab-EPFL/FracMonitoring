function [ind] = aic_pick(x,o)
% 
% comment by Dmitry Loginov 23/09/2019: more or less the same we use
% 

%    written by Dong Liu 17/09/2019
%    Modified from  Kalkan, E.'s Matlab codes.
%    P-PHASE ARRIVAL TIME PICKER BASED ON AKAIKE INFORMATION CRITERION
%
%    Computes P-phase arrival time in windowed digital single-component
%    acceleration or broadband velocity record without requiring threshold
%    settings. Returns P-phase arrival time index.
%
%    An abbreviated form of the Akaike Information Criterion is used
%    since we are interested in the global minimum. See paper by
%    Maeda, N. (1985). A method for reading and checking phase times in
%    autoprocessing system of seismic wave data, ZisinJishin 38, 365?379.
%
%   Syntax:
%          ind = aic_pick(x,'whole') consider whole seismogram to be picked.
%
%          ind = aic_pick(x,'to_peak') consider part of seismogram from 1
%          to peak value of x.
%
%   Input:
%          x = raw data in single-column format
%
%   Output:
%          ind = Index of pick in x; pick is the min(AIC(n)) + 1
%
%          AIC(n) = k*log(var(x[1,k])) + (n-k-1)*log(var(x[k+1,n])) where k
%          goes from 1 to length(x) 

x = x - median(x); % remove median of x and window
switch o
    case {0,'1','to_peak'}
        ind_peak = find(abs(x) == max(abs(x)));
        xnew = x(1:ind_peak);
    otherwise
        xnew = x;
end

junk = aicval(xnew);

if junk ~= 0;
    ind = find(junk == min(junk)) + 1; % pick is one more than divide point
else
    ind = 0;
end

% show the plot of the AIC index and together the signal
figure
plot(1:size(junk,2),(junk)')
xamp=max(junk)-min(junk);
hold on 
plot(1:size(xnew,1),xnew/max(xnew)*xamp+mean(junk))
hold on
line([ind ind], [min(junk) max(junk)],'Color', 'black','LineStyle', '--')

% show the plot of the AIC index and together the signal, sccaled onto [-1,1]
% xnewtest=(junk-mean(junk))/xamp;
% figure
% plot(1:size(junk,2),xnewtest')
% hold on 
% plot(1:size(xnew,1),xnew/max(xnew))
% hold on
% line([ind ind], [min(xnewtest) max(xnewtest)],'Color', 'black','LineStyle', '--')


return

function [a] = aicval(x)
if ~isempty(x)
    n = length(x);
    for i=1:n-1;
        %compute variance in first part
        s1 = var(x(1:i));
        if s1 <= 0;
            s1 = 0;
        else
            s1=log(s1);
        end
        %compute variance in second part
        s2 = var(x(i+1:n));
        if s2 <= 0;
            s2 = 0;
        else
            s2=log(s2);
        end
        a(i) = i*(s1) + (n-i+1)*(s2);
    end
else
    a = 0;
end
return