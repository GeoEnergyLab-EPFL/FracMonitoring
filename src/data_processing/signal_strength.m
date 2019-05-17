function [normAll, h_SNR]  = signal_strength(dataset,endnoise,varargin)
% function computing basic signal strength for all pairs and extracting the
% "diagonal" corresponding to on-axis pairs, with optional matrix plot.
% varargin can be endnoise, the length of excitation noise at the beginning
% of the recorded signals

% open figure from passed handle if it exists
if ~isempty(varargin)
    for i_arg = 1:nargin-1
        if isgraphics(varargin{1})
            figure(varargin{1});
            hold on
        else
            h_SNR = figure;
        end
    end
else
    h_SNR = figure;
end


% compute L2-norm strength
normAll = squeeze(vecnorm(dataset(endnoise:end,:,:),2,1));

% plot
figure(h_SNR)
imagesc(1:size(dataset,3),1:size(dataset,2),normAll)
axis square
caxis([0 1]*2.5)
colorbar
colormap('jet')
xlabel('Receiver #')
ylabel('Source #')

% % EXPERIMENTAL: SNR ESTIMATION
% % compute L2-norm of the early noise section
% noiseAll = squeeze(vecnorm(dataset(endnoise:endflat,:,:),2,1));
% SNRAll = normAll./noiseAll;
% 
% % plot
% figure
% imagesc(1:size(dataset,3),1:size(dataset,2),SNRAll)
% axis square
% caxis([0 1]*4)
% colormap('jet')
% xlabel('Receiver #')
% ylabel('Source #')

end