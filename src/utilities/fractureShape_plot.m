function [fig_handle] = fractureShape_plot(seq_i,m,Solid,SRPairs,myBlock,myTransducers,myPlattens,AcqTime,varargin)
% Dong Liu -- 10/10/2019
% plot the ellipse fracture with diffracted points and SR pairs and the
% corresponding opening if there is any
% the last two arguments decide whether to plot the fracture opening
% argument 1: the fig handle
% argument 2: the path for opening sequence number
% argument 3: the path for opening
% the last two arguments can be easily adjusted when the format of the opening
% changes

ell = Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
res = diffractionForward(Solid,SRPairs,ell);% give one the shortest time needed for diffraction

narg = length(varargin);
if narg>=1
    if ~isempty(varargin)
                if isgraphics(varargin{1})
                    fig_handle = figure(varargin{1});
                else
                    fig_handle = figure;
                end
     else
                fig_handle = figure;
     end
     hold on
end

% plot the corresponding width profile
plotblockwithplattens(myBlock,myPlattens,fig_handle)
plotdiffrays(SRPairs,res(:,2),res(:,3),res(:,4),fig_handle);% one should plot the trays with the corresponding diffracted points
plotEllipse(ell,fig_handle,'b.-');
hold on
plot3(res(:,2),res(:,3),res(:,4),'.g','MarkerSize',30);
title(['\fontsize{40}Seq ' num2str(seq_i) ': ' datestr(AcqTime(seq_i))])

% fracture opening plot
if narg>=3
    if ~isempty(varargin) && ~isempty(varargin{2}) && ~isempty(varargin{3})
        SRtrsn = SourceReceiverPairs(myTransducers,myPlattens,[1:16;1:16]');
        trsn_x = SRtrsn.XS_XR(1:end,1); % x-coordinate
        trsn_y = SRtrsn.XS_XR(1:end,2); % y-coordinate
        trsn_z = (SRtrsn.XS_XR(1:end,3)+SRtrsn.XS_XR(1:end,6))/2; % z-coordinate
        % load the global sequence number
        [seq_list] = importdata(varargin{2},'\t');
        % load the fracture opening
        [width_profile] = importdata(varargin{3},'\t');
        [~,idx] = ismember(seq_i,seq_list');% we can do this since the picked-sequence is the global sequence
        width_picked = (width_profile(idx,1:end))'/max(max(width_profile))*0.02;
        % 0.02 here is just for plotting the amplitude of opening
        hold on
        widthplot1 = fill3([trsn_x(1); trsn_x(1:8); trsn_x(8)],[trsn_y(1);trsn_y(1:8); trsn_y(8)],[trsn_z(1);trsn_z(1:8)+width_picked(1:8);trsn_z(8)],'r');
        alpha(widthplot1,0.2)
        set(widthplot1,'EdgeColor','none')
        hold on
        widthplot2 = fill3([trsn_x(9); trsn_x(9:16); trsn_x(16)],[trsn_y(9);trsn_y(9:16); trsn_y(16)],[trsn_z(9);trsn_z(9:16)+width_picked(9:16);trsn_z(16)],'r');
        alpha(widthplot2,0.2)
        set(widthplot2,'EdgeColor','none')
    else
        disp('Wrong path for fracture opening, there will be no fracture opening plotted')
    end
end

end