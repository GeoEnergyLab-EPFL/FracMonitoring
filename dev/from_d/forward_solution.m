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

debug = 'yes'; % 'yes' for debug or 'no' for simple execution

%% define material properties for TI material

% FORWARD_SOLID
[TIalpha, TIbeta, TIvec, rho, c11, c33, c55, c66, c13, npts, t, UP, Psi] = forward_solid(debug)


%% define source and receiver positions

% FORWARD_ELLIPSE
[xc, yc, zc, r1, r2, Ealpha, Ebeta, Egamma, Elli0, EulerRot, Elli1, Ev1, Ev2, En]...
    = forward_ellipse(debug, npts, t)

%% define source and receiver positions

% FORWARD_SR
[ptop, pbottom, pnorth, psouth, peast, pwest, ...
    dcross, dhoriz, dvert1, dvert2, ddiag1, ddiag2,...
    Xsources, Ysources, Zsources, Sources, Xreceivers, Yreceivers, Zreceivers, Receivers,...
    SStmp, RRtmp, SS, RR]...
    = forward_SR(debug, Elli1)

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

% mark source and receiver with black dor %%% mm?
switch debug
    case 'yes'
        figure(f_main)
        plot3(Sources(1,ss),Sources(2,ss),Sources(3,ss),'.k')
        plot3(Receivers(1,rr),Receivers(2,rr),Receivers(3,rr),'.k')
end

% line containing both source and receiver
Lv = Sources(:,ss)-Receivers(:,rr); %%% ?

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

% find points at min dist+20% %%% 20?
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
% find group velocity for that angle %%% what is UP?
UpRD = Up(PsiRDidx);
% compute travel-time
ttRD = dRD./UpRD;

% find point on ellipse with min total travel-time
[TT, Emintt] = min(ttSD+ttRD);

results(:,qq) = [TT; Emintt];

% plot on figure
switch debug
    case 'yes'
        figure%(f_main)
        plot3(Elli1(1,Inear(Emintt)),Elli1(2,Inear(Emintt)),Elli1(3,Inear(Emintt)),'oc')
end
end

% get execution time
toc