clearvars
home

% create block object
myBlock = SampleBlock([250 250 250]);

% create two platten objects
PlattenA = Platten('A','B',myBlock);
PlattenC = Platten('C','W',myBlock);

%figA = plattenplot2D(PlattenA);
%figC = plattenplot2D(PlattenC);

fig3D = plattenplot3D(PlattenA,myBlock);
plattenplot3D(PlattenC,myBlock,fig3D);

% parameters for tranducer object
transd_serial = [1 2 3 4 5 6 7 8 10 11];
transd_type = {'P' 'P' 'P' 'P' 'S' 'P' 'S' 'S' 'P' 'P'};
transd_channel = {'S1' 'R2' 'S3' 'R4' 'S5' 'R6' 'S7' 'R8' 'R20' 'S21'};
transd_locid = [8 9 10 11 12 13 14 15 3 6];
transd_platten = {'A' 'A' 'A' 'A' 'A' 'A' 'A' 'A' 'C' 'C'};
% create transducer object
myTransd = Transducers(transd_serial,transd_type,transd_channel,transd_platten,transd_locid);

% get transducer coordinates
xyz = myTransd.calc_global_coord({PlattenA,PlattenC});
plot3(xyz(:,1),xyz(:,2),xyz(:,3),'b.')