
npts = 10; % how many points to sample the original curve. 
t = linspace(0,2*pi,npts);
%z = linspace(-1,1,npts);
a=0.3;b=1;
xyz_toapp = [cos(t).*a; sin(t).*b; t*0.]; % create a XYZ matrix of samples point on the curve to spline - make sure the first and last point are the same to have close contour


spell=cscvn(xyz_toapp(:,[1:end 1])); % call matlab spline function 

% checkin that it's ok
hold on
npts = 220;
t = linspace(0,2*pi,npts);
xyz = [cos(t).*a; sin(t).*b; t*0.];

plot3(xyz(1,:),xyz(2,:),xyz(3,:),'ro','LineWidth',2); axis equal;
fnplt(spell,'k',2)
hold off


% function to get the XYZ of a point along the curve knowing its
% curvilinear axis
% the curvilinear axis bounds are the first and last entry of the breaks
% field of the spline  object (here it is spell)
ppval(spell,[0,1])

% gett the "central" point of the curve

mean(ppval(spell,spell.breaks)')




