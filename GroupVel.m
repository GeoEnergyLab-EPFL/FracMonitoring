function [Up, Usv, Ush, Psi] = GroupVel(c11,c33,c13,c55,c66,rho,theta)
%GROUPVEL computes the group velocities of a TI medium for discrete angles
%   GROUPVEL returns the group velocities and group andle of a TI medium as
%   a function of the elastic constants C11, C33, C13, C55, C66, the
%   density RHO, for a vector of phase angles theta
%
% Thomas Blum, Geo-Energy Lab, EPFL, November 2018

% Equations describing propagation in TI medium from Tsvankin (2001) book

% symbolic derivations of the P-wave group velocity
% define symbolic variables
syms thetaSym
tmp1Sym = (c11+c55)*sin(thetaSym).^2 + (c33+c55)*cos(thetaSym).^2;
tmp2Sym = (c11-c55)*sin(thetaSym).^2 - (c33-c55)*cos(thetaSym).^2;
tmp3Sym = 4*(c13+c55)^2*sin(thetaSym).^2.*cos(thetaSym).^2;

% compute phase velocity
VpSym = sqrt((tmp1Sym*1E9+sqrt((tmp2Sym*1E9).^2+tmp3Sym*1E18))/(2*rho));

% and derivative
dVpSym = diff(VpSym,thetaSym);
% compute group velocity (Eq 1.70 in Tsvankin)
UpSym = VpSym.*sqrt(1+(1./VpSym.*dVpSym).^2);

% estimate for a range of phase angles
Up = double(subs(UpSym,thetaSym,theta));

% shear velocities (not yet computed)
Ush = NaN;
Usv = NaN;

% group angle (Eq. 1.71 in Tsvankin)
num = double(subs(dVpSym/VpSym,thetaSym,theta));
denom = sin(theta).*cos(theta).*(1-tan(theta).*num);
tanPsi = tan(theta).*(1+num./denom);
Psi = atan(tanPsi);

% fix NaN (from division by zero)
Psi(isnan(Psi)) = 0;

% custom-made unwrap
dPsi = diff(Psi);
% fix positive jumps
Ipos = find(dPsi>=0.9*pi);
for ii = Ipos
    Psi(ii+1:end) = Psi(ii+1:end)-pi;
end
% fix negative jumps
Ineg = find(dPsi<=-0.9*pi);
for ii = Ineg
    Psi(ii+1:end) = Psi(ii+1:end)+pi;
end

end