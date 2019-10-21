function [ind] = pai_kpicker(x,M,o)
%    written by Dong Liu 17/09/2019
%    P-PHASE ARRIVAL TIME PICKER BASED ON phase arrival
%    identification-kurtosis
%
%    Computes P-phase arrival time with moving windows. 
%    Returns P-phase arrival time index.
%
%   Input:
%          x = raw data in single-column format
%          M = the window size, should be smaller than the signal length
%   Output:
%          ind = Index of pick in x; pick is the min(AIC(n)) + 1
%

x = x - median(x); % remove median of x and window
switch o
    case {0,'1','to_peak'}
        ind_peak = find(abs(x) == max(abs(x)));
        xnew = x(1:ind_peak);
    otherwise
        xnew = x;
end
junk=zeros(size(xnew,1),1);
for k=M:size(xnew,1)
    junk(k,1) = kfunc(xnew,k,M);
end
junk(1:M-1)=min(junk);
if junk ~= 0
    ind = find(junk == max(junk)); % corresponds to the maximum value
else
    ind = 0;
end

% show the plot of the AIC index and together with the signal
figure
plot(1:size(junk,1),junk)
xamp=max(junk)-min(junk);
hold on
plot(1:size(xnew,1),xnew/max(xnew)*xamp+mean(junk))
hold on
line([ind ind], [min(junk) max(junk)],'Color', 'black','LineStyle', '--')

return

function [a] = kfunc(x,k,M)
if ~isempty(x)
    xnew=x(k-M+1:k);
    mk=mean(xnew);
    aup=0;
    adown=0;
    for i=k-M+1:k
        aup=aup+(x(i)-mk)^4;
        adown=adown+(x(i)-mk)^2;
    end
    a=(M-1)*aup/adown-3;
else
    a = 0;
end
return