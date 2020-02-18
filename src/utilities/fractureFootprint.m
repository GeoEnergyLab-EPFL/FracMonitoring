function [fig_handle] = fractureFootprint(m,seqlist,SRPairs_all,ray_type_all, myBlock,varargin)
% Dong Liu -- 28/11/2019
% plot the ellipse fracture with diffracted points 
% argument 1: the fig handle
% argument 2: colorstyle for ellipse or radial
% argument 3: plotstyle for ellipse or radial

global m_ind;
global z_c;
global Solid;

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
plotstyle='-';
colorstyle=[1 0 0];
if narg>=2 && ~isempty(varargin) && ~isempty(varargin{2})
    colorstyle=varargin{2};
    if narg>=3 && ~isempty(varargin{2})
        plotstyle=varargin{3};
    end
end

nseq=length(seqlist);
% 2D fracture plot function
% input: sequence list, plot style
for i = 1:nseq
        m_i = m(seqlist(i),:);
        hold on
        rectangle('Position',[0 0 myBlock.L_N myBlock.L_E]);
        hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%         switch length(m_i)
%             case 8 % Ellipse
%                 ell = Ellipse(m_i(1),m_i(2),m_i(3:5),m_i(6),m_i(7),m_i(8));
%                 [Elli_pts]=PointsOnEllipse(ell,720); % discretizing the ellipse with 120 pts
%                 
%             case 7 % Ellipse with fixed z_c
%                 ell = Ellipse(m_i(1),m_i(2),[m_i(3:4),z_c],m_i(5),m_i(6),m_i(7));
%                 [Elli_pts]=PointsOnEllipse(ell,720); % discretizing the ellipse with 120 pts
%                 %plot(m_i(3),m_i(4),'r.','MarkerSize',10);
%             case 6
%                 ell = Radial(m_i(1),m_i(2:4),m_i(5),m_i(6));
%                 [Elli_pts]=PointsOnRadial(ell,720); % discretizing the ellipse with 120 pts
%                 %plot(m_i(2),m_i(3),'r.','MarkerSize',10);
%             case 5
%                 ell = Radial(m_i(1),[m_i(2:3),z_c],m_i(4),m_i(5));
%                 [Elli_pts]=PointsOnRadial(ell,720); % discretizing the ellipse with 120 pts
%                 %plot(m_i(2),m_i(3),'r.','MarkerSize',10);
%             otherwise
%                 disp('Please check your input m-vector')
%                 return;
%         end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
       switch m_ind
            case 1 % Ellipse
                ell = Ellipse(m_i(1),m_i(2),m_i(3:5),m_i(6),m_i(7),m_i(8));
                [Elli_pts]=PointsOnEllipse(ell,720); % discretizing the ellipse with 120 pts
            case 2
                ell = Radial(m_i(1),m_i(2:4),m_i(5),m_i(6));
                [Elli_pts]=PointsOnRadial(ell,720); % discretizing the ellipse with 120 pts
            case 3 % planar Ellipse
                ell = Ellipse(m_i(1),m_i(2),m_i(3:5),m_i(6),0,0);
                [Elli_pts]=PointsOnEllipse(ell,720); % discretizing the ellipse with 120 pts
            case 4 % planar Radial
                ell = Radial(m_i(1),m_i(2:4),0,0);
                [Elli_pts]=PointsOnRadial(ell,720); % discretizing the ellipse with 120 pts                
            case 5 % Ellipse with fixed z_c
                ell = Ellipse(m_i(1),m_i(2),[m_i(3:4),z_c],m_i(5),m_i(6),m_i(7));
                [Elli_pts]=PointsOnEllipse(ell,720); % discretizing the ellipse with 120 pts
            case 6
                ell = Radial(m_i(1),[m_i(2:3),z_c],m_i(4),m_i(5));
                [Elli_pts]=PointsOnRadial(ell,720); % discretizing the ellipse with 120 pts
            otherwise
                disp('Please check your input m-vector')
                return;
        end
        
        hold on
        plot(Elli_pts(:,1),Elli_pts(:,2),plotstyle,'color',colorstyle,'LineWidth',2.5);
        
        hold on
        xlim([0, myBlock.L_N]);
        ylim([0, myBlock.L_E]);
        xlabel('Easting (m)');
        ylabel('Northing (m)');

        res = diffractionForward(Solid,SRPairs_all(seqlist(i)),ell,ray_type_all{seqlist(i)});
        plot(res(:,2),res(:,3),'.','color',[252 241 155]/255.,'MarkerSize',4);
        
end
hold on
if (m_ind==1) || (m_ind==3) || (m_ind==5) % the case of the ellipse
    plot(m(seqlist,3),m(seqlist,4),'.-','color',colorstyle,'MarkerSize',4);
else
    plot(m(seqlist,2),m(seqlist,3),'.-','color',colorstyle,'MarkerSize',4);
end

end