function [results]=diffractionForward(Solid,SRPairs,EllipseObj)

% description:
% compute the diffraction travel time for a list of sources receivers
% pair knowing a diffractor curve.
%
% inputs:
% Solid : solid properties object either Isotropic or anisotropic
% SRPairs: SR pairs map object containing all necessary info (notably the
% diffracted wave type)
% EllipseObj : ellipse object for the diffractor or a radial object for the
% diffractor
%
% Dt_diffracted: corresponding diffracted arrival time and location XYZ
%  of the diffracting point
%
% Brice Lecampion - 2019

% isotropy only !
%
% note for extension to TI - be careful as in the block coordinate system
% the medium is not VTI ! 

x_source=SRPairs.XS_XR(:,1:3);
x_rec=SRPairs.XS_XR(:,4:6);

% En=EllipseObj.Normal;

classname=class(EllipseObj);
switch classname 
    case 'Ellipse'
    Elli1=PointsOnEllipse(EllipseObj,720);
    case 'Radial'
    Elli1=PointsOnRadial(EllipseObj,720);
    otherwise
        disp('the input must be an ellipse or raidal object');
end

if isa(Solid,'IsotropicSolid')
    
    vP=Solid.Vp;
    results=zeros(SRPairs.n_pairs,4);
    
    % vectorize everything -> 300 much faster !
    dsd=pdist2(x_source,Elli1);
    ddr=pdist2(x_rec,Elli1);
    dt=(dsd+ddr)/vP;
    [res, ks]=min(dt,[],2);
    results(:,1)=res;
    results(:,2:4)=Elli1(ks,:);
    
else
    disp('Error only Isotopric solid for now');
    return
end
end