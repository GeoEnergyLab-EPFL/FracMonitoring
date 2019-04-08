function [normAll, diagonal, h_SNR]  = signal_strength(dataset,varargin)
% function computing basic signal strength for all pairs and extracting the
% "diagonal" corresponding to on-axis pairs, with optional matrix plot.
% varargin can be endnoise, the l

% open figure from passed handle if it exists
if ~isempty(varargin)
    for i_arg = 1:nargin-1
    if isnumeric(varargin{1})
        figure(varargin{1});
        hold on
    else
        fig_handle = figure;
    end
else
    fig_handle = figure;
end

end