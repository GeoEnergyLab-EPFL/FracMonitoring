function fig_handle = plot_cuboid(lengths,origin,varargin)
% PLOTCUBOID - Display a 3D-cuboid in a new or existing figure
%
%   PLOTCUBE(LENGTHS,ORIGIN) displays a 3D-cuboid in a new graphical object
%   with the following properties
%   * LENGTHS : 3-elements vector that defines the length of cuboid edges
%   * ORIGIN: 3-elements vector that defines the start point of the cube
%   (from lower left bottom corner)
%
%   PLOTCUBE(LENGTHS,ORIGIN,FIGHANDLE) displays the cuboid in the existing
%   graphical object with the handle FIG_HANDLE

% create vertices for unit cube
tmp1 = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
% stretch to proper size
tmp2 = tmp1.*lengths(:)';
% shift to proper origin
xyzVertices = tmp2+origin(:)';
% create face sequence for patch
faces = [1 2 3 4 1 5 6 2 6 7 3 7 8 4 8 5];

% open figure from passed handle if it exists
if ~isempty(varargin)
    if isgraphics(varargin{1})
        figure(varargin{1});
        hold on
    else
        fig_handle = figure;
    end
else
    fig_handle = figure;
end
% plot the cuboid
patch('Faces',faces,'Vertices',xyzVertices,'EdgeColor','black','FaceColor',...
    'none','LineWidth',2)
view(3)
end