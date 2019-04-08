
function [xc, yc, zc, r1, r2, Ealpha, Ebeta, Egamma, Elli0, EulerRot, Elli1, Ev1, Ev2, En]...
    = forward_ellipse(debug, npts, t)

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
En = cross(Ev1,Ev2); %%% ?

end
