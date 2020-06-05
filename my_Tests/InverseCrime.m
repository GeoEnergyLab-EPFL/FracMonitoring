% forward test script 

%% cleanup first, set global parameters
close all
clearvars
home

% data storage location
datastor = 'gel-nas1'; % 'local' if dataset copied to local drive, 'gel-nas1', or 'enacdrives'

%% choose dataset and load acquisition times
switch datastor
    case 'gel-nas1'
        % data path on gel-nas1
        datapath = pathbyarchitecture('gel-nas1');
    case 'enacdrives'
        datapath = pathbyarchitecture('enac1files');
    case 'local'
        [~, username] = system('whoami');
        datapath = ['/home/' username(1:end-1) '/data/'];
end

% 2019 acquisitions
datayear = 19;
% % injection on slate
datamonth = 03;
dataday = 01;
starttime = '121330';
% test on gabbro
%datamonth = 03;
%dataday = 22;
%starttime = '145034';

% data folder name from experiment date
datafold = [num2str(datayear,'%02d') '-' num2str(datamonth,'%02d') '-' ...
    num2str(dataday,'%02d') '/'];

% extract header info from JSON file
fjson = [datapath datafold num2str(starttime) '.json'];
[jsonhdr,myTransducers,myPlattens,myBlock] = load_header(fjson);


 
%%  Prepare for diffraction forward - select some S-R pairs

mySR1=TwoPlattenPairs(myTransducers,myPlattens,'B','W') ; 
mySR2=TwoPlattenPairs(myTransducers,myPlattens,'B','E') ; 
mySR3=TwoPlattenPairs(myTransducers,myPlattens,'B','N') ; 
mySR4=TwoPlattenPairs(myTransducers,myPlattens,'B','S') ; 

the_comb_map=CombineMaps(mySR1,mySR2);

the_comb_map2=CombineMaps(mySR3,mySR4);
my_map=union(the_comb_map,the_comb_map2,'rows');

my_map=[27 7; 
    23 8];

SRdiff=SourceReceiverPairs(myTransducers,myPlattens,my_map(:,:));
fig_b=plotblockwithplattens(myBlock,myPlattens)

%
fig_handle=plotdirectrays(SRdiff,fig_b);

%%
%----- 
marble = IsotropicSolid(2700,60.*1e9,0.25); % vel in m/s
%%%
mtrue=[90/1000;40/1000; 0.25/2;.25/2;.25/2;0*pi/180;45*pi/180;0.*pi/180];
elltrue = Ellipse(mtrue(1),mtrue(2),mtrue(3:5),0*pi/180,45*pi/180,0.*pi/180);
fig_handle=plotEllipse(elltrue,fig_handle);

%% ----
 

%%
tic
[results]=diffractionForward(marble,SRdiff,elltrue);
toc
%% plot things
fig_b=plotblockwithplattens(myBlock,myPlattens)

%
fig_handle=plotdirectrays(SRdiff,fig_b);
fig_handle=plotEllipse(elltrue,fig_handle);

hold on
plot3(results(:,2),results(:,3),results(:,4),'.g','MarkerSize',30);

%%

dclean=results(:,1);

% noise
vv=0.03*var(dclean)
gm = gmdistribution(0.,vv);

noise=random(gm,length(dclean));

plot(dclean,'.-'); hold on;
plot(dclean+noise,'.r'); hold on;

%%

global d  SRPairs 
d=dclean+noise;
SRPairs=SRdiff;

global Solid  Cdinvdiag;
Solid=marble;
Cdinvdiag = 1/(1*vv).*(0.*dclean+1); % data inverse of variance

global prior

mp= [ .025;.025; .125;.125;.125; 0.;45*pi/180.;0. ];
sig_p=[.125;.125; 0.01;0.01;0.01;0.5;0.05;0.5];
prior=GaussianPrior(mp,sig_p);

%%
test = minusPosteriorPDF(1.*mp)


%% using a differential evolution algorithm
VTR=0.1;
D=length(mp);
XVmin=0.*mp';
XVmax=[.250;.250;.250;.250;.250;90;90;180]';
NP = 10*D;
itermax = 1000; 
F = 0.8; 
CR = 0.5; 
strategy = 7;
refresh = 10; 
[bestmem,bestval,nfeval] = devec3(@minusPosteriorPDF,VTR,D,XVmin,XVmax,[],NP,itermax,F,CR,strategy,refresh);


%% 
m=bestmem';
ell=Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
res= diffractionForward(Solid,SRPairs,ell);

plot(d); hold on;
plot(res(:,1));


%%
fig_b=plotblockwithplattens(myBlock,myPlattens)

fig_handle=plotdirectrays(SRdiff,fig_b);

fig_handle=plotEllipse(elltrue,fig_handle);


fig_handle=plotEllipse(ell,fig_handle,'b.-');

hold on
plot3(res(:,2),res(:,3),res(:,4),'.g','MarkerSize',30);


%% MCMC

mstart=m;
Cstart=diag(1./prior.invCpdiag);
Smax=20000*length(m);
fid_log = fopen('test.txt','w');

[accept, Xf, PXf, cf ,k_o, stab]=MarkovChainMonteCarlo(fid_log, mstart,Cstart,Smax,1,0,@minusPosteriorPDF);


%%
acceptrate=mean(accept) % should be between 0.2-0.4


% resampling step
ns=1;
q=length(Xf(1,:)); % problem dimension


tmp_ndx = find(accept);

k_s=[k_o:ns:length(PXf)];
 
 
mLogP=PXf(k_s);
Xsp=Xf(k_s,:);


% -LOG POST histogram
handler_histpost=figure('Name','Histogram of - Log Posterior');
hist(mLogP,10);  % it should look like a lognormal pdf
title('Histogram of - Log Posterior');
legend('It should look like a log Normal pdf' );

%%

m=mean(Xsp);

m=mean(Xf(20000:30000,:));

ell=Ellipse(m(1),m(2),m(3:5),m(6),m(7),m(8));
res= diffractionForward(Solid,SRPairs,ell);

plot(d); hold on;
plot(res(:,1));


fig_b=plotblockwithplattens(myBlock,myPlattens)

%
fig_handle=plotdirectrays(SRdiff,fig_b);
fig_handle=plotEllipse(elltrue,fig_handle);
fig_handle=plotEllipse(ell,fig_handle,'b');

hold on
plot3(res(:,2),res(:,3),res(:,4),'.g','MarkerSize',30);

%%

%plot(Xsp(:,1),Xsp(:,2),'.')

plot(Xsp(:,8),Xsp(:,7),'.')


%%
%%[MPost,Ctilde_all]=ProcessingMCMC(accept,PX,X,k_o,stab,mcmc_res_logfile)

% [MPost,Ctilde_all]=ProcessingMCMC(accept,PXf,Xf,k_o,stab,@minusPosteriorPDF,[]);
 