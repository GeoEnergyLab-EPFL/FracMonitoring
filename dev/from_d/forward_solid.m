
function [TIalpha, TIbeta, TIvec, rho, c11, c33, c55, c66, c13, npts, t, Up, Psi]...
    = forward_solid(debug)

% define axis of TI anisotropy (vector)
TIalpha = 0*pi/180; % angle between TI axis and x-axis
TIbeta = 90*pi/180; % angle between TI axis and horizontal plane
TIvec = [cos(TIalpha)*cos(TIbeta); sin(TIalpha)*cos(TIbeta); sin(TIbeta)]; %%% why?

% marble
rho_marble = 2700;  % density marble (from litterature)
%  elastic constants in GPa
c13_marble = 0; % lambda in isotropy case %%% whatis lambda and mu?
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


% discretize range of angles
npts = 720; % nb of points %%% what is it?
t = linspace(0,2*pi,npts+1);
t = t(1:end-1);
% estimate for a range of phase angles
[Up, ~, ~, Psi] = GroupVel(c11,c33,c13,c55,c66,rho,t);

% check with plot %%% what is the check?
switch debug
    case 'yes'
        f_main = figure
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
end
