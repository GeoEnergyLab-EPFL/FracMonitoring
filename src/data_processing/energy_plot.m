function [F,fig_handle]  = energy_plot(myEnergy,myBlock,myPlattens,myTransducers,DiffractionRecord,varargin)
% Dong Liu -- 18/09/2019
% function creating plots on energy change compared with the reference signal for
% all transmisson transducer pairs on the top and bottom platten or for the
% source or receiver and its corresponding receivers or sources on the
% opposite platten
% optional argument 1: fig_handle
% optional argument 2: sequence number and Time (not coded)
% optional argument 2: color index for fracture front or fluid front
% optional argument 3: criteria to triger the bi-color color change
% DiffractionRecord can be empty or must have at least four properties
% .seqnb,  .mDE, and .acqT and .m_ind are least required.


if ~isempty(varargin) 
    if isgraphics(varargin{1})
        fig_handle = figure(varargin{1});
    else
        fig_handle = figure;
    end
else
    fig_handle = figure;
end

if ~isempty(varargin) && ~isempty(varargin{2})
    colortype=varargin{2};
else
    colortype='r';
end

if ~isempty(varargin) && ~isempty(varargin{3})
    yield=varargin{3};
else
    yield=[];
end

if ~isempty(myEnergy) && ~isempty(myEnergy.seq)
    SRmap=myEnergy.SRmap;
    seq=myEnergy.seq;
    energy=myEnergy.energy;
else
    return
end

hold on

if SRmap(1:end,1) == SRmap(1:end,2)
    [F,fig_handle] = evolution_transmit(seq,SRmap,energy,myBlock,myPlattens,myTransducers,DiffractionRecord,fig_handle,colortype,yield);
else
    [F,fig_handle] = evolution_single(seq,SRmap,energy,myBlock,myPlattens,myTransducers,DiffractionRecord,fig_handle,colortype,yield);
end

end

function [F,fig_handle] = evolution_transmit(seq,SRmap,energy,myBlock,myPlattens,myTransducers,DiffractionRecord,fig,colortype,yield)
sidemarker='T'; % we plot the top platten by default 
fig_handle=plotside2Dwithplattens(myBlock,myPlattens,sidemarker,fig);
hold on
mkrsize=8;
xyzTransd = calc_global_coord(myTransducers,myPlattens);

for i=1:size(energy,1)
    plotside2Dwithplattens(myBlock,myPlattens,sidemarker,fig_handle);
    
    % add the sequence number and the acquisition time
    if isempty(DiffractionRecord)
        if length(seq)==size(energy_all,1)
            title(['Seq ' num2str(seq(i))]);
        end
    else 
        % there is sequences recorded
        % plot the ellipse
        idx=0;
        for i_d=1:size(DiffractionRecord.mDE,1)
            if idx==0
            if DiffractionRecord.seqnb(i_d)==seq(i)
            	idx=i_d;
            end
            end
        end
        disp(idx)
        if idx>0
            m=DiffractionRecord.mDE(idx,:);
            % we need to ajust the model 
            m_ind=DiffractionRecord.m_ind;
            switch m_ind
                case 1
                    ell=Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
                    [Elli_pts]=PointsOnEllipse(ell,720); % discretizing the ellipse with 120 pts
                case 2
                    ell=Radial(m(1),m(2:4),m(5),m(6));
                    [Elli_pts]=PointsOnRadial(ell,720); % discretizing the ellipse with 120 pts
                case 3
                    ell=Ellipse(m(1),m(2),m(3:5),m(6),0,0);
                    [Elli_pts]=PointsOnEllipse(ell,720); % discretizing the ellipse with 120 pts
                case 4
                    ell=Radial(m(1),m(2:4),0,0);
                    [Elli_pts]=PointsOnRadial(ell,720); % discretizing the ellipse with 120 pts
                case 5
                    disp('Fixed z-corordinate ellipse is not coded yet');
                case 6
                    disp('Fixed z-corordinate circle is not coded yet');
                otherwise
                    disp('Please check the m_ind input in the Diffraction Record file');
            end
            
            plot(Elli_pts(:,1),Elli_pts(:,2),[colortype '-'],'LineWidth',2.5);
            hold on
            %title(['Seq ' num2str(DiffractionRecord.seqnb(idx)) ' ' DiffractionRecord.acqT(idx)])
            hold on
        end
    end  
    
    % plot all the corresponding transducers and their energy evolution
    % ratio
    for ii=1:size(SRmap)
    % set tge switch trigger value
        if ~isempty(yield)
                % set color value    
            if energy(i,ii)>=yield  %7 for G01 fracture opening  %0.8 for M05 transmission signal
                clr=[1 1 1];
                %clr=[253 231 37]/255.;
            else
                clr=[253 231 37]/255.;
                %clr=[1 1 1];
            end
            hold on
            plot(xyzTransd(SRmap(ii,1),1),xyzTransd(SRmap(ii,1),2),...
        'o','MarkerSize',mkrsize,'MarkerFaceColor',clr,'MarkerEdgeColor','k');
        else
% get coordinate from the S or R transducers, they share the same
% coordinates on x and y
            hold on;
            plot(xyzTransd(SRmap(ii,1),1),xyzTransd(SRmap(ii,1),2),...
         'o','MarkerSize',mkrsize,'MarkerFaceColor',...
         [1 1 1].*heaviside(1-energy(i,ii))*abs(energy(i,ii)-1)/max(abs(energy(1:end,ii)-1)),...
         'MarkerEdgeColor','k');
                
        end
    end
    
    F(i)=getframe(fig_handle);
    pause(0.1);
    if i<size(energy,1)
        clf(gcf);
    end
end

end


function [F,fig_handle] = evolution_single(seq,SRmap,energy_all,myBlock,myPlattens,myTransducers,DiffractionRecord,fig,colortype,yield)
% now only working for the top and bottom platten
% we locate the transducer position with the number and the type
if size(SRmap,1)==1
    number=SRmap(1,1);
    type='S';
elseif SRmap(1,1)==SRmap(2,1)
    number=SRmap(1,1);
    type='S';
elseif SRmap(1,2)==SRmap(2,2)
    number=SRmap(1,2);
    type='R';
else
    disp('Wrong input for SRmap, only single Platten groups accepted')
end

id_transducer = find((myTransducers.channel+1==number).*(myTransducers.type==type));

for p=1:length(myPlattens)
    if myTransducers.platten(id_transducer)==myPlattens(p).id
        oposide=myPlattens(p).face;
    end
end
switch oposide
    case 'N'
        sidemarker='S';
    case 'S'
        sidemarker='N';
    case 'E'
        sidemarker='W';
    case 'W'
        sidemarker='E';
    case 'T'
        sidemarker='B';
    case 'B'
        sidemarker='T';
end 

if type=='S'
    k=2;
else
    k=1;
end

fig_handle = plotside2Dwithplattens(myBlock,myPlattens,sidemarker,fig);
hold on
% plot the all the receivers if source, all the sources if receiver
% we build the SR Pairs
% we calculate the energy relative change
% we set these changes as the color map
mkrsize=18;
switch sidemarker
    case 'N'
        i=1;
        j=3;
    case 'S'
        i=1;
        j=3;
    case 'E'
        i=2;
        j=3;
    case 'W'
        i=2;
        j=3;
    case 'T'
        i=1;
        j=2;
    case 'B'
        i=1;
        j=2;
    otherwise
        disp('Wrong input for side indicator')

end

xyzTransd = calc_global_coord(myTransducers,myPlattens);

% if it is source k=2; if it is receiver k=1;
for s=1:size(energy_all,1)
    plotside2Dwithplattens(myBlock,myPlattens,sidemarker,fig_handle);
    
    % set the title and plot the diffracted ellipse if there is any
    if isempty(DiffractionRecord)
        if length(seq)==size(energy_all,1)
            title(['Seq ' num2str(seq(s))]);
        end
    else 
        % there is sequences recorded
        % plot the ellipse
        idx=0;
        for i_d=1:size(DiffractionRecord.mDE,1)
            if DiffractionRecord.seqnb(i_d)==seq(s)
            	idx=i_d;
            end
        end
        disp(idx)
        if idx>0
            m=DiffractionRecord.mDE(idx,:);
            % we need to ajust the model 
            m_ind=DiffractionRecord.m_ind;
            switch m_ind
                case 1
                    ell=Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
                    [Elli_pts]=PointsOnEllipse(ell,720); % discretizing the ellipse with 120 pts
                case 2
                    ell=Radial(m(1),m(2:4),m(5),m(6));
                    [Elli_pts]=PointsOnRadial(ell,720); % discretizing the ellipse with 120 pts
                case 3
                    ell=Ellipse(m(1),m(2),m(3:5),m(6),0,0);
                    [Elli_pts]=PointsOnEllipse(ell,720); % discretizing the ellipse with 120 pts
                case 4
                    ell=Radial(m(1),m(2:4),0,0);
                    [Elli_pts]=PointsOnRadial(ell,720); % discretizing the ellipse with 120 pts
                case 5
                    disp('Fixed z-corordinate ellipse is not coded yet');
                case 6
                    disp('Fixed z-corordinate circle is not coded yet');
                otherwise
                    disp('Please check the m_ind input in the Diffraction Record file');
            end
            
            plot(Elli_pts(:,i),Elli_pts(:,j),[colortype '-'],'LineWidth',3);
            hold on
            title(['Seq ' num2str(DiffractionRecord.seqnb(idx)) ' ' DiffractionRecord.acqT(idx)])
            hold on
        end
    end

    for ii=1:size(SRmap,1)
    % get coordinates from the SRmap
        idx = find(myTransducers.channel==SRmap(ii,k)-1);
        if ~isempty(yield)
            if energy(s,ii)>=yield  %7 for G01 fracture opening  %0.8 for M05 transmission signal
                clr=[1 1 1];
            else
                clr=[1 0 0];
            end
            hold on
            plot(xyzTransd(idx,i),xyzTransd(idx,j),...
        'o','MarkerSize',mkrsize,'MarkerFaceColor',clr,'MarkerEdgeColor','k');
        else
             hold on
            plot(xyzTransd(idx,i),xyzTransd(idx,j),...
        'o','MarkerSize',mkrsize,'MarkerFaceColor',...
        [1 1 1].*heaviside(1-energy_all(s,ii))*abs(energy_all(s,ii)-1)/max(abs(energy_all(1:end,ii)-1)),...
        'MarkerEdgeColor','k');
            
        end
       
    end
    hold on
    F(s)=getframe(fig_handle);
    pause(0.1);
    if s<size(energy_all,1)
        clf(fig_handle);
    end
end

end