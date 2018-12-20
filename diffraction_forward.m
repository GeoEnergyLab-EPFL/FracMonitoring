%% Script for tip diffraction forward problem
% assumes elliptical fracture shape
% all positions and distances in mm, reference is the center of the block
% frame of reference
% x axis - EW
% y axis - NS
% z axis - vertical

close all
clearvars
home

debug = 'no'; % 'yes' for debug or 'no' for simple execution

%% define material properties for TI material
% as a first step we look at P-waves only

% define axis of TI anisotropy (vector)
TIalpha = 0*pi/180; % angle between TI axis and x-axis
TIbeta = 90*pi/180; % angle between TI axis and horizontal plane
TIvec = [cos(TIalpha)*cos(TIbeta); sin(TIalpha)*cos(TIbeta); sin(TIbeta)];

% marble
rho_marble = 2700;  % density marble (from litterature)
%  elastic constants in GPa
c13_marble = 0; % lambda in isotropy case
c66_marble = 0; % mu in isotrpy case
c11_marble = 0; % lambda+2mu
c33_marble = 0; % lambda+2mu
c55_marble = 0; % mu

% carmen slate
rho_carmen = 2770;
c11_carmen = 100;
c33_carmen = 50;
c55_carmen = 20;
c66_carmen = 40;
c13_carmen = 40;

% choose the right solid properties
solid = 'carmen';
rho = eval(['rho_' solid]);
c11 = eval(['c11_' solid]);
c33 = eval(['c33_' solid]);
c55 = eval(['c55_' solid]);
c66 = eval(['c66_' solid]);
c13 = eval(['c13_' solid]);

%% compute TI group velocities as a function of group angle for a set of angles
% discretize range of angles
npts = 720; % nb of points
t = linspace(0,2*pi,npts+1);
t = t(1:end-1);
% estimate for a range of phase angles
[Up, ~, ~, Psi] = GroupVel(c11,c33,c13,c55,c66,rho,t);

% check with plot
switch debug
    case 'yes'
        figure
        % plot group velocity
        yyaxis left
        plot(t*180/pi,Up,'*')
        xlabel('Phase angle (\circ)')
        ylabel('Group velocity (m/s)')
        % plot angles
        hold on
        yyaxis right
        plot(t*180/pi,Psi*180/pi,'or')
        ylabel('Group angle (\circ)')
        axis tight
end


%% define source and receiver positions
% platen position wrt center of block
ptop = 125;
pbottom = -125;
pnorth = 125;
psouth = -125;
peast = 125;
pwest = -125;

% spacing on platen
dcross = 18;    % top platen cross spacing
dhoriz = 50;    % side platen horiz spacing
dvert1 = 20;    % side platen first vert spacing
dvert2 = 30;    % side platen farther vert spacing
ddiag1 = 34;    % top platen first diag spacing
ddiag2 = 28;    % top platen additional diag spacing

% position of source transducers (at center)
Xsources = [zeros(1,8), [4 3 2 1]*dcross, -[1 2 3 4]*dcross,...
    ones(1,6)*peast, -dhoriz , 0, dhoriz, dhoriz, -dhoriz, 0,...
    [ddiag1+ddiag2 ddiag1 -ddiag1 -ddiag1-ddiag2]/sqrt(2)];
Ysources = [[4 3 2 1]*dcross, -[1 2 3 4]*dcross, zeros(1,8),...
     dhoriz , 0, -dhoriz, -dhoriz, 0, dhoriz, ones(1,6)*pnorth,...
     [ddiag1+ddiag2 ddiag1 ddiag1 ddiag1+ddiag2]/sqrt(2)];
Zsources = [ones(1,16)*ptop, ones(1,3)*-dvert1, ones(1,3)*dvert1,...
    ones(1,3)*-(dvert1+dvert2), ones(1,3)*(dvert1+dvert2),...
    ones(1,4)*ptop];
% combine all three coordinates
Sources = [Xsources; Ysources; Zsources];

% position of receiver transducers (at center)
Xreceivers = [Xsources(1:16), ones(1,6)*pwest, Xsources(23:end)];
Yreceivers = [Ysources(1:22), ones(1,6)*psouth, Ysources(29:end)];
Zreceivers = [ones(1,16)*pbottom, Zsources(17:28), ones(1,4)*pbottom];
% combine all three coordinates
Receivers = [Xreceivers; Yreceivers; Zreceivers];


%% ellipse parameters
% center position
xc = 0;
yc = 0;
zc = 0;

% long and short axes
r1 = 75;
r2 = 50;

% Euler angles for ellipse
Ealpha = 30*pi/180;
Ebeta = 0*pi/180;
Egamma = 0*pi/180;

%% compute ellipse
% ellipse before rotation using discretized angles
Elli0 = [xc+r1*cos(t); yc+r2*sin(t); zc+0*t];

% Euler rotation matrix zxz
EulerRot = [cos(Ealpha)*cos(Egamma)-cos(Ebeta)*sin(Ealpha)*sin(Egamma),...
    -cos(Ealpha)*sin(Egamma)-cos(Ebeta)*cos(Egamma)*sin(Ealpha),...
    sin(Ealpha)*sin(Ebeta);...
    cos(Egamma)*sin(Ealpha)+cos(Ealpha)*cos(Ebeta)*sin(Egamma),...
    cos(Ealpha)*cos(Ebeta)*cos(Egamma)-sin(Ealpha)*sin(Egamma),...
    -cos(Ealpha)*sin(Ebeta);...
    sin(Ebeta)*sin(Egamma), cos(Egamma)*sin(Ebeta), cos(Ebeta)];

% rotate
Elli1 = EulerRot*Elli0;

% extract vectors for plane containing ellipse
Ev1 = Elli1(:,1);
Ev2 = Elli1(:,round(npts/2));
% normal vector
En = cross(Ev1,Ev2);


%% plot sources and receivers
switch debug
    case 'yes'
        f_main = figure;
        % plot source locations
        plot3(Xsources,Ysources,Zsources,'or')
        % plot receiver locations
        hold on
        plot3(Xreceivers,Yreceivers,Zreceivers,'ob')
        % add ellipse
        plot3(Elli1(1,:),Elli1(2,:),Elli1(3,:),'k')
        % fix figure
        axis equal
        xlabel('Easting')
        ylabel('Northing')
        zlabel('Depth')
end


%% compute for series of source-receiver pairs instead
% sources on top
SStmp = [5*ones(1,6) 6*ones(1,6) 7*ones(1,6) 8*ones(1,6) 13*ones(1,6)...
    14*ones(1,6) 15*ones(1,6) 16*ones(1,6)];
% receivers on the side
RRtmp = [[23 24 25 26 27 28], [23 24 25 26 27 28], [23 24 25 26 27 28],...
    [23 24 25 26 27 28], [17 18 19 20 21 22], [17 18 19 20 21 22],...
    [17 18 19 20 21 22], [17 18 19 20 21 22]];
% and vice-versa, combined
SS = [SStmp RRtmp];
RR = [RRtmp SStmp];

%% forward propagation: first locate candidates for diffraction point
% choose source and receiver pair
%ss = 7;
%rr = 27;

% measure execution time
tic
% allocate result matrix
results = zeros(2,length(SS));
% loop on source-receiver pairs instead
for qq = 1:length(SS) % possible to parralellize with parfor
    ss = SS(qq);
    rr = RR(qq);

% mark source and receiver with black dor
switch debug
    case 'yes'
        figure(f_main)
        plot3(Sources(1,ss),Sources(2,ss),Sources(3,ss),'.k')
        plot3(Receivers(1,rr),Receivers(2,rr),Receivers(3,rr),'.k')
end

% line containing both source and receiver
Lv = Sources(:,ss)-Receivers(:,rr);

% intersection of line and ellipse plane
t = dot(En,[xc; yc; zc]-Sources(:,ss))/dot(En,Lv);
Intersec = Sources(:,ss)+t*Lv;

% add intersection to plot
switch debug
    case 'yes'
        figure(f_main)
        plot3(Intersec(1),Intersec(2),Intersec(3),'g*')
end

% search for closest point on ellipse
% compute distances
dElli = dists(Intersec,Elli1);

% find points at min dist+20%
Inear = find(dElli<=min(dElli)*1.2);
% add to plot
switch debug
    case 'yes'
        figure(f_main)
        plot3(Elli1(1,Inear),Elli1(2,Inear),Elli1(3,Inear),'ok')
end

% forward propagation: find propagation times for these candidates
% compute the distances between source, diffractor and receiver
dSD = dists(Sources(:,ss),Elli1(:,Inear));
dRD = dists(Receivers(:,rr),Elli1(:,Inear));
% compute source-diffractor, and diffractor-receiver direction vectors
directSD = directions(Sources(:,ss),Elli1(:,Inear));
directRD = directions(Receivers(:,rr),Elli1(:,Inear));

% add graphical check to main figure
switch debug
    case 'yes'
        % vector scaling
        vscal = 4;
        % add vectors
        figure(f_main)
        quiver3(Elli1(1,Inear),Elli1(2,Inear),Elli1(3,Inear),...
            directSD(1,:),directSD(2,:),directSD(3,:),vscal,'r')
        quiver3(Elli1(1,Inear),Elli1(2,Inear),Elli1(3,Inear),...
            directRD(1,:),directRD(2,:),directRD(3,:),vscal,'b')
end

% source-diffractor (SD) travel-time
% find group angles
PsiSD = atan2(vecnorm(bsxfun(@cross,TIvec,directSD)),dot(repmat(TIvec,1,...
    length(Inear)),directSD));
% find indices for group angles
PsiSDidx = zeros(size(PsiSD));
for ii = 1:length(Inear)
    [~, tmpidx] = min(abs(Psi-PsiSD(ii)));
    PsiSDidx(ii) = tmpidx;
end
% find group velocities for these angles
UpSD = Up(PsiSDidx);
% compute travel-time
ttSD = dSD./UpSD;

% receiver-diffractor (RD) travel-time
% find group angle
PsiRD = atan2(vecnorm(bsxfun(@cross,TIvec,directRD)),dot(repmat(TIvec,1,...
    length(Inear)),directRD));
% find indices for group angles
PsiRDidx = zeros(size(PsiRD));
for ii = 1:length(Inear)
    [~, tmpidx] = min(abs(Psi-PsiRD(ii)));
    PsiRDidx(ii) = tmpidx;
end
% find group velocity for that angle
UpRD = Up(PsiRDidx);
% compute travel-time
ttRD = dRD./UpRD;

% find point on ellipse with min total travel-time
[TT, Emintt] = min(ttSD+ttRD);

results(:,qq) = [TT; Emintt];

% plot on figure
switch debug
    case 'yes'
        figure(f_main)
        plot3(Elli1(1,Inear(Emintt)),Elli1(2,Inear(Emintt)),Elli1(3,Inear(Emintt)),'oc')
end
end

% get execution time
toc