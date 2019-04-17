function [C, lags] = xcorr1point5(A,B)

% "1.5D" crosscorrelation between A, a 2D array and a column vector B, of
% same length as the colums of A.

if ~isequal(size(A,1),length(B))||~isvector(B)
    fprintf('\nError: uncompatible input sizes!\n');
    return
end

% lengths
np = length(B); % number of points in time series
nt = size(A,2); % number of traces

[C, lags] = xcorr([A, B]);
C = reshape(C,2*np-1,nt+1,nt+1);
C = flipud(C(:,1:nt,nt+1));

end
