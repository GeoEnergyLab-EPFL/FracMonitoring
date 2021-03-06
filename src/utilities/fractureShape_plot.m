function [fig_handle] = fractureShape_plot(m,Solid,SRPairs,ray_type, myBlock,myTransducers,myPlattens,varargin)
% Dong Liu -- 10/10/2019
% plot the ellipse fracture with diffracted points and SR pairs and the
% corresponding opening if there is any
% the last two arguments decide whether to plot the fracture opening
% argument 1: the fig handle
% argument 2: plotstyle for ellipse or radial
% argument 3: the path for opening sequence number
% argument 4: the path for opening
% argument 5: the global sequence number you want to plot the fracture opening  
% the last three arguments can be easily adjusted when the format of the opening
% changes

global m_ind;
global z_c;

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
% set the default plotstyle for the fracture geometry
plotstyle='b.-';
if narg>=2 && ~isempty(varargin) && ~isempty(varargin{2})
    plotstyle=varargin{2};
end

plotblockwithplattens(myBlock,myPlattens,fig_handle);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% switch length(m)
%     case 8 % Ellipse
%         ell = Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
%         plotEllipse(ell,fig_handle,plotstyle);
%     case 7 % Ellipse with fixed z_c
%         ell = Ellipse(m(1),m(2),[m(3:4),z_c],m(5),m(6),m(7));
%         plotEllipse(ell,fig_handle,plotstyle);
%     case 6
%         ell = Radial(m(1),m(2:4),m(5),m(6));
%         plotRadial(ell,fig_handle,plotstyle);
%     case 5
%         ell = Radial(m(1),[m(2:3),z_c],m(4),m(5));
%         plotRadial(ell,fig_handle,plotstyle);
%     otherwise
%         disp('Please check your input m-vector')
%         return;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch m_ind
    case 1 % Ellipse
        ell = Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
        plotEllipse(ell,fig_handle,plotstyle);
    case 2
        ell = Radial(m(1),m(2:4),m(5),m(6));
        plotRadial(ell,fig_handle,plotstyle);
    case 3 % Ellipse
        ell = Ellipse(m(1),m(2),m(3:5),m(6),0,0);
        plotEllipse(ell,fig_handle,plotstyle);
    case 4
        ell = Radial(m(1),m(2:4),0,0);
        plotRadial(ell,fig_handle,plotstyle);    
    case 5 % Ellipse with fixed z_c
        ell = Ellipse(m(1),m(2),[m(3:4),z_c],m(5),m(6),m(7));
        plotEllipse(ell,fig_handle,plotstyle);
    case 6
        ell = Radial(m(1),[m(2:3),z_c],m(4),m(5));
        plotRadial(ell,fig_handle,plotstyle);
    otherwise
        disp('Please check your input m-vector')
        return;
end


% we set this by default for the PdP
% ray_type=ones(size(SRPairs.SRmap,1),1);
res = diffractionForward(Solid,SRPairs,ell,ray_type);% give one the shortest time needed for diffraction
plotdiffrays(SRPairs,res(:,2),res(:,3),res(:,4),fig_handle);% one should plot the trays with the corresponding diffracted points
hold on
plot3(res(:,2),res(:,3),res(:,4),'.g','MarkerSize',30);


% fracture opening plot
if narg>=4
    if ~isempty(varargin) && ~isempty(varargin{3}) && ~isempty(varargin{4}) && ~isempty(varargin{5})
        SRtrsn = SourceReceiverPairs(myTransducers,myPlattens,[1:16;1:16]');
        trsn_x = SRtrsn.XS_XR(1:end,1); % x-coordinate
        trsn_y = SRtrsn.XS_XR(1:end,2); % y-coordinate
        trsn_z = (SRtrsn.XS_XR(1:end,3)+SRtrsn.XS_XR(1:end,6))/2; % z-coordinate
        % load the global sequence number
        [seq_list] = importdata(varargin{3},'\t');
        % load the fracture opening
        [width_profile] = importdata(varargin{4},'\t');
        seq_i = varargin{5}; % the sequence number at which you want to plot the opening
        [~,idx] = ismember(seq_i,seq_list');
        if idx>0
            width_picked = (width_profile(idx,1:end))'/max(max(width_profile))*0.02;
        % 0.02 here is just for plotting the amplitude of opening
        else
            width_picked = zeros(16,1);
        end
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