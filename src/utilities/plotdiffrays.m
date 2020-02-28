function fig_handle = plotdiffrays(obj,res_x,res_y,res_z,varargin)
        % Dong Liu -- 09/09/2019
        % res_x: vector of x-coordinates of the diffracted points
        % res_y: vector of y-coordinates of the diffracted points
        % res_z: vector of z-coordinates of the diffracted points
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
            
            for i=1:obj.n_pairs
                plot3([obj.XS_XR(i,1) res_x(i)]',[obj.XS_XR(i,2) res_y(i)]',[obj.XS_XR(i,3) res_z(i)]','.-.','color',[0.5,0.5,0.5],'linewidth',2)
                hold on;
                plot3([res_x(i) obj.XS_XR(i,4)]',[res_y(i) obj.XS_XR(i, 5)]',[res_z(i) obj.XS_XR(i,6)]','.-.','color',[0.5,0.5,0.5],'linewidth',2)
                hold on;
            end
            plot3(obj.XS_XR(:,1),obj.XS_XR(:,2),obj.XS_XR(:,3),'r.','MarkerSize',18);
            hold on
            plot3(obj.XS_XR(:,4),obj.XS_XR(:,5),obj.XS_XR(:,6),'b.','MarkerSize',18);
            
            %
            axis equal
end